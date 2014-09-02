function varargout = plot(varargin)
%PLOT  Plots the feasible region of a set of constraints
%
% p = plot(C,x,c,n,options)
%
% Note that only convex sets are allowed, or union of convex sets
% represented using binary variables (either defined explictly or
% introduced by YALMIP when modelling, e.g., mixed integer linear
% programming representable operators)
%
% C:  Constraint object
% x:  Plot variables [At most three variables]
% c:  color [double] ([r g b] format) or char from 'rymcgbk'
% n:  #vertices [double ] (default 100 in 2D and 300 otherwise)
% options: options structure from sdpsettings

% Get the onstraints
if nargin<1
    return
end

F = varargin{1};

if length(F)==0
    return;
end

if nargin < 5
    opts = sdpsettings('verbose',0);
else
    opts = varargin{5};
    if isempty(opts)
        opts = sdpsettings('verbose',0);
    end
end
opts.verbose = max(opts.verbose-1,0);

if any(is(F,'uncertain'))
    F = robustify(F,[],opts);
end

if nargin < 3
    color=['rymcgbk']';
else
    color = varargin{3};
    if isa(color,'sdpvar')
        error('The variables should be specified in the second argument.');
    end
    color = color(:)';
    if isempty(color)
        color=['rymcgbk']';
    end
end

% Plot onto this projection (at most in 3D)
if nargin < 2
    x = [];
else
    x = varargin{2};
    if ~isempty(x)
        if ~isa(x,'sdpvar')
            error('Second argument should be empty or an SDPVAR');
        end
        x = x(:);
        x = x(1:min(3,length(x)));
    end    
end

if isempty(F)
    return
end

if any(is(F,'sos'))
    % Use image representation, feels safer (I'm not sure about the logic
    % in the code at the moment. Can the dualization mess up something
    % otherwise...)
    if ~(opts.sos.model == 1)
        opts.sos.model = 2;
    end
    % Assume the variables we are plotting are the parametric
    F = compilesos(F,[],opts,x);
end

% Create a model in YALMIPs low level format
% All we change later is the cost vector
%sol = solvesdp(F,sum(x),opts);
[model,recoverdata,diagnostic,internalmodel] = export(F,[],opts,[],[],0);
if isempty(internalmodel) | (~isempty(diagnostic) & diagnostic.problem)
    error('Could not create model. Can you actually solve problems with this model?')
end
internalmodel.options.saveduals = 0;
internalmodel.getsolvertime = 0;
internalmodel.options.dimacs = 0;

% Try to find a suitable set to plot
if isempty(x)
    if isempty(internalmodel.extended_variables) & isempty(internalmodel.aux_variables)
        x = depends(F);
        x = x(1:min(3,length(x)));
        localindex = 1;
        localindex = find(ismember(recoverdata.used_variables,x));
    else
        % not extended variables
        candidates = setdiff(1:length(internalmodel.c),[ internalmodel.aux_variables(:)' internalmodel.extended_variables(:)']);
        % Not nonlinear variables
        candidates = candidates(find(internalmodel.variabletype(candidates)==0));
        % Settle with this guess
        localindex = candidates(1:min(3,length(candidates)));
        x = localindex;
    end
else
    localindex = [];
    for i = 1:length(x)
        localindex = [localindex find(ismember(recoverdata.used_variables,getvariables(x(i))))];
    end
end

if nargin < 4
    if length(x)==3
        n = 300;
    else
        n = 100;
    end
else
    n = varargin{4};
    if isempty(n)
        if length(x)==3
            n = 300;
        else
            n = 100;
        end
    end
    if ~isa(n,'double')
        error('4th argument should be an integer>0');
    end
end

if ~isempty(internalmodel.integer_variables)
    error('PLOT can currently not display sets involving integer variables');
end

if isempty(internalmodel.binary_variables)
    [x_opt{1},errorstatus] = generateBoundary(internalmodel,x,n,localindex);
else
    if strcmp(internalmodel.solver.tag,'BNB')
        internalmodel.solver = internalmodel.solver.lower;
    end
    nBin = length(internalmodel.binary_variables);
    p = extractLP(internalmodel);
    p = extractOnly(p,internalmodel.binary_variables);
    
    internalmodel.F_struc = [zeros(nBin,size(internalmodel.F_struc,2));internalmodel.F_struc];
    I = eye(nBin);
    internalmodel.F_struc(1:nBin,1+internalmodel.binary_variables) = I;
    internalmodel.K.f = internalmodel.K.f + length(internalmodel.binary_variables);
    errorstatus = 1;
    x_opt = {};
    errorstatus = zeros(1,2^nBin);
    for i = 0:2^nBin-1;
        comb = dec2decbin(i,nBin);
        if checkfeasiblefast(p,comb(:),1e-6)
            internalmodel.F_struc(1:nBin,1) = -comb(:);
            [x_temp,wrong] = generateBoundary(internalmodel,x,n,localindex);
            if ~wrong
                errorstatus(i+1) = 0;
                x_opt{end+1} = x_temp;
            end
        else
            errorstatus(i+1)=0;
        end
    end
end

if all(errorstatus)
    if nargout==0
        plot(0);
    end
elseif nargout == 0
    for i = 1:length(x_opt)
        try
            plotSet(x_opt{i},color(1+rem(i-1,size(color,1)),:),opts);
        catch
        end
    end
end

if nargout > 0
    varargout{1} = x_opt;
end


function [xout,errorstatus] = solvefordirection(c,internalmodel,uv)
internalmodel.c = 0*internalmodel.c;
internalmodel.c(uv) = c;
sol  = feval(internalmodel.solver.call,internalmodel);
xout = sol.Primal;
xout = xout(uv(:));
errorstatus = sol.problem;


function p = plotSet(x_opt,color,options)
if size(x_opt,1)==2
    p = patch(x_opt(1,:),x_opt(2,:),color);
    set(p,'LineStyle',options.plot.wirestyle);   
    set(p,'LineWidth',options.plot.linewidth);
    set(p,'EdgeColor',options.plot.edgecolor);
    set(p,'Facealpha',options.plot.shade);    
else
    try
        K = convhulln(x_opt');
        p = patch('Vertices', x_opt','Faces',K,'FaceColor', color);
    catch
         p = fill3(x_opt(1,:),x_opt(2,:),x_opt(3,:),1);
    end
    set(p,'LineStyle',options.plot.wirestyle);   
    set(p,'LineWidth',options.plot.linewidth);
    set(p,'EdgeColor',options.plot.edgecolor);
    set(p,'Facealpha',options.plot.shade);     
    lighting gouraud;
    view(3);
    camlight('headlight','infinite');
    camlight('headlight','infinite');
    camlight('right','local');
    camlight('left','local');   
end




function [x_opt,errorstatus] = generateBoundary(internalmodel,x,n,localindex);

x_opt = [];
phi = [];
errorstatus = 0;
waitbar_created = 0;
t0 = clock;
waitbar_starts_at = 2;
lastdraw = clock;
try % Try to ensure that we close h
    if length(x)==2
        mu = 0.5;
    else
        mu=1;
    end
    n_ = n;
    n = ceil(mu*n);
   % h = waitbar(0,['Please wait, solving ' num2str(n_) ' problems using ' internalmodel.solver.tag]);
    angles = (0:(n))*2*pi/n;
    if length(x)==2
        c = [cos(angles);sin(angles)];
    else
        c = randn(3,n);
    end
    i=1;
    while i<=n & errorstatus ~=1
        [xi,errorstatus] = solvefordirection(c(:,i),internalmodel,localindex(:));
        if errorstatus == 2
            disp('Discovered unbounded direction. You should add bounds on variables')            
        elseif errorstatus == 12
            [xi,errorstatus] = solvefordirection(0*c(:,i),internalmodel,localindex(:));
            if errorstatus == 0
                errorstatus = 2;
                disp('Discovered unbounded direction. You should add bounds on variables')
            end
        end                                        
        x_opt = [x_opt xi];
        if ~waitbar_created
            if etime(clock,t0)>waitbar_starts_at;
             h = waitbar(0,['Please wait, solving ' num2str(n_) ' problems using ' internalmodel.solver.tag]);
             waitbar_created = 1;
            end
        end
        if waitbar_created & etime(clock,lastdraw)>1/10
            waitbar(i/n_,h)
            lastdraw = clock;
        end
        i=i+1;
    end
    
    if errorstatus==0 & length(x)==2
        % Close the set
        x_opt = [x_opt x_opt(:,1)];
        
        % Add points adaptively
        pick = 1;
        n = floor((1-mu)*n_);
        for i = 1:1:n
            for j= 1:(size(x_opt,2)-1)
                d = x_opt(:,j)-x_opt(:,j+1);
                distance(j,1) = d'*d;
            end
            [dist,pos]=sort(-distance);
            % Select insertion point
            phii=(angles(pos(pick))+angles(pos(pick)+1))/2;
            xi = solvefordirection([cos(phii);sin(phii)],internalmodel,localindex);
            d1=xi-x_opt(:,pos(pick));
            d2=xi-x_opt(:,pos(pick)+1);
            if d1'*d1<1e-3 | d2'*d2<1e-3
                pick = pick+1;
            else
                angles = [angles(1:pos(pick)) phii  angles((pos(pick))+1:end)];
                x_opt = [x_opt(:,1:pos(pick)) xi  x_opt(:,(pos(pick))+1:end)];
            end
            if ~waitbar_created
            if etime(clock,t0)>waitbar_starts_at;
                h = waitbar(0,['Please wait, solving ' num2str(n_) ' problems using ' internalmodel.solver.tag]);
                waitbar_created = 1;
            end
            end
            if waitbar_created
            waitbar((ceil(n_*mu)+i)/n_,h);
            end
        end
    end
    if waitbar_created        
    close(h);
    end
catch
    if waitbar_created
    close(h);
    end
end

function pLP = extractLP(p);
pLP = p;
pLP.F_struc = pLP.F_struc(1:p.K.f+p.K.l,:);
pLP.K.q = 0;
pLP.K.s = 0;

function pRed = extractOnly(p,these);
pRed = p;
p.F_struc(:,1+these) = 0;
removeEQ = find(any(p.F_struc(1:pRed.K.f,2:end),2));
removeLP = find(any(p.F_struc(1+pRed.K.f:end,2:end),2));
pRed.F_struc(pRed.K.f+removeLP,:)=[];
pRed.F_struc(removeEQ,:)=[];
pRed.K.f = pRed.K.f - length(removeEQ);
pRed.K.l = pRed.K.l - length(removeLP);
pRed.F_struc = pRed.F_struc(:,[1 1+these]);
pRed.lb = pRed.lb(these);
pRed.ub = pRed.ub(these);