function [g,geq,dg,dgeq,xevaled] = fmincon_con(x,model,xevaled)

global latest_xevaled
global latest_x_xevaled
global sdpLayer

% Early bail for linear problems
g    = [];
geq  = [];
dg   = [];
dgeq = [];

% Nothing nonlinear, and no cones
if model.linearconstraints && ~any(model.K.q) && ~any(model.K.e) && ~any(model.K.s)
    xevaled = [];
    return
end

% Compute all nonlinear variables
if nargin<3
    if isequal(x,latest_x_xevaled)
        xevaled = latest_xevaled;
    else
        xevaled = zeros(1,length(model.c));
        xevaled(model.linearindicies) = x;
        xevaled = apply_recursive_evaluation(model,xevaled);
        latest_x_xevaled = x;
        latest_xevaled = xevaled;
    end
end

% All simple scalar inequalities
if model.nonlinearinequalities
    g = full(model.Anonlinineq*xevaled(:)-model.bnonlinineq);    
end

% Append SOCP cones
if any(model.K.q)
    top = startofSOCPCone(model.K);
    z0 = model.F_struc*[1;xevaled];
    for i = 1:length(model.K.q)
        z = z0(top:top+model.K.q(i)-1);
        g = [g;-(z(1) - sqrt(z(2:end)'*z(2:end)))];
        top = top + model.K.q(i);
    end
end

% Append EXP cones 
if any(model.K.e)
    top = startofEXPCone(model.K);
    % x3-x2*exp(x1/x2)>=0
    for i = 1:model.K.e        
        X = model.F_struc(top:top+2,:)*[1;xevaled];        
        if X(2)==0
            % Eeew, nasty closure let's hope for the best
            r = X(3);
        else
            r = X(3) - X(2)*exp(X(1)/X(2));
        end
        g = [g;-r];
        top = top + 3;          
    end
end

% Append POW cones 
if any(model.K.p)
    top = startofPOWCone(model.K);
    % x1^a x2^(1-a)>=norm(x(:end))
    for i = 1:length(model.K.p)
        n = model.K.p(i);
        alpha = model.F_struc(top+n-1,1);
        X = model.F_struc(top:top+n-2,:)*[1;xevaled];       
        r = X(1)^alpha*X(2)^(1-alpha)-norm(X(3:end));
        g = [g;-r];
        top = top + n;          
    end
end

% Append SDP cones (eigenvalues)
if any(model.K.s)
    top = startofSDPCone(model.K);
    for i = 1:length(model.K.s)
        n = model.K.s(i);
        X = model.F_struc(top:top+n^2-1,:)*[1;xevaled];
        X = full(reshape(X,n,n));
        [d,v] = eig(X); 
        v = diag(v);
        if strcmpi(model.options.slayer.algorithm,'convex') || isequal(model.options.slayer.algorithm,1)
            v = cumsum(v);
        end
        % These will reordered later
        g = [g;-v(1:min(n,sdpLayer.n))];
        top = top + n^2;          
    end
end

% Create all nonlinear equalities
if model.nonlinearequalities
    geq = full(model.Anonlineq*xevaled(:)-model.bnonlineq);
end

% Now Jacobians...
if nargout == 2 || ~model.derivative_available
    return
elseif isempty(model.evalMap) && (model.nonlinearinequalities==0) && (model.nonlinearequalities==0) && (model.nonlinearcones==0) && anyCones(model.K)
    % Linear SOCP
    dg_SOCP = computeSOCPDeriv(model,xevaled);
    % Linear EXP
    dg_EXP = computeEXPDeriv(model,xevaled);
    % Linear POW
    dg_POW = computePOWDeriv(model,xevaled);
    % Linear SDP
    dg_SDP = computeSDPDeriv(model,xevaled,1);
    dg = [dg_SOCP;dg_EXP;dg_POW;dg_SDP];     
    dg = dg(:,model.linearindicies);      
    g = reorderEigenvalueG(g,model); 
    
elseif isempty(model.evalMap) & (model.nonlinearinequalities || model.nonlinearequalities || model.nonlinearcones) 
    % Nonlinear terms in constraints so jacobian stuff needed
    % however, there are only monomials, which is exploited
    newdxx = computeMonomialVariableJacobian(model,x);
        
    if ~isempty(model.Anonlineq)      
        dgeq = model.Anonlineq*newdxx; 
    end
    if ~isempty(model.Anonlinineq)
        dg = model.Anonlinineq*newdxx;                
    end    
    if any(model.K.q)
        dg = [dg;computeSOCPDeriv(model,xevaled,newdxx)];
    end  
    if any(model.K.e)
        dg = [dg;computeEXPDeriv(model,xevaled,newdxx)];
    end 
    if any(model.K.p)
        dg = [dg;computePOWDeriv(model,xevaled,newdxx)];
    end      
    if any(model.K.s)
        dg_SDP = computeSDPDeriv(model,xevaled,newdxx);               
        dg = [dg;dg_SDP]; 
        g = reorderEigenvalueG(g,model);       
    end
else  
    % Completely general case with monomials and functions
    requested = model.fastdiff.requested;
    dx = apply_recursive_differentiation(model,xevaled,requested,model.Crecursivederivativeprecompute);        
    if ~isempty(model.Anonlineq)
        dgeq = [model.Anonlineq*dx];  
    end
    if ~isempty(model.Anonlinineq)
        dg = [model.Anonlinineq*dx];
    end        
    if any(model.K.q)
        dg = [dg;computeSOCPDeriv(model,xevaled,dx)];
    end    
    if any(model.K.e)
        dg = [dg;computeEXPDeriv(model,xevaled,dx)];
    end 
    if any(model.K.p)
        dg = [dg;computePOWDeriv(model,xevaled,dx)];
    end     
    if any(model.K.s)
        dg_SDP = computeSDPDeriv(model,xevaled,dx);                       
        dg = [dg;dg_SDP];
        g = reorderEigenvalueG(g,model);       
    end
end

% For performance reasons, we work with transpose
if model.nonlinearequalities
    dgeq = dgeq';
end
if model.nonlinearinequalities || anyCones(model.K)
    dg = dg';   
end
% hold on
% x = sdpvar(2,1);
% plot([g + dg'*(x-xevaled(1:2))<=0,-3<=x<=3],x,'y',[],sdpsettings('plot.shade',.1))
% plot(xevaled(1),xevaled(2),'+b')
% drawnow
% axis([-2 2 -2 2])
% dg
% xevaled
% 1;

function conederiv = computeSOCPDeriv(model,z,dzdx)
conederiv = [];
z = sparse(z(:));
if any(model.K.q)
    top = startofSOCPCone(model.K);
    for i = 1:length(model.K.q)
        d = model.F_struc(top,1);
        b = model.F_struc(top+1:top+model.K.q(i)-1,1);
        cA = model.F_struc(top:top+model.K.q(i)-1,2:end);
        c = cA(1,:)';
        A = cA(2:end,:);
        % -c'*x - d + ||Ax+b||<=0
        if nargin == 2
            % No inner derivative
            if length(model.linearindicies)==length(model.c)
                e = A*z + b;
                smoothed = sqrt(10^-10 + e'*e);
                temp = -c+(A'*(b + A*z))/smoothed;
                conederiv = [conederiv temp];
            else
                A = A(:,model.linearindicies);
                c = c(model.linearindicies);                
                e = A*z(model.linearindicies) + b;
                smoothed = sqrt(10^-10 + e'*e);
                temp = (-c'+(A'*(b + A*z(model.linearindicies)))'/smoothed);
                conederiv = [conederiv temp'];
            end
        else
            % -c'*x - d + ||Ax+b||>=0
            e = A*z + b;
            
            smoothed = sqrt(10^-10 + e'*e);
            temp = (-c'+(A'*(b + A*z))'/smoothed);
            conederiv = [conederiv (temp*dzdx)'];
        end
        top = top + model.K.q(i);
    end
end
conederiv = conederiv';

function conederiv = computeEXPDeriv(model,z,dzdx)
conederiv = [];
z = sparse(z(:));
if any(model.K.e)
    top = startofEXPCone(model.K);
    for i = 1:model.K.e
        c1 = model.F_struc(top,2:end);
        c2 = model.F_struc(top+1,2:end);
        c3 = model.F_struc(top+1,2:end);
        x = model.F_struc(top:top+2,:)*[1;z];
        x1 = x(1);x2 = x(2);x3 = x(3);
        % x3 - x2exp(x1/x2) >= 0
        % d/dx = [-exp(x1) -(exp(x1/x2)-x1/x2*exp(x1/x2)) 1]       
        if nargin == 2
            c1=c1(model.linearindicies);
            c2=c2(model.linearindicies);
            c3=c3(model.linearindicies);
            temp = c3-(c2*exp(x1/x2)+exp(x1/x2)*(c1*x2-c2*x1)/x2);
            conederiv = [conederiv -temp'];
        else
            temp = c3-(c2*exp(x1/x2)+exp(x1/x2)*(c1*x2-c2*x1)/x2);
            conederiv = [conederiv -(temp*dzdx)'];
        end
        top = top + 3;
    end
end
conederiv = conederiv';

function conederiv = computePOWDeriv(model,z,dzdx)
conederiv = [];
z = sparse(z(:));
if any(model.K.p)
    top = startofPOWCone(model.K);
    for i = 1:length(model.K.p)
        n = model.K.p(i);
        c1 = model.F_struc(top,2:end);
        c2 = model.F_struc(top+1,2:end);
        c3 = model.F_struc(top+2:top+n-2,2:end);
        alpha = model.F_struc(top+n-1,1);
        x = model.F_struc(top:top+n-2,:)*[1;z];
        x1 = x(1);x2 = x(2);x3 = x(3:end);
        % x1^alpha*x2^(1-alpha) - norm(x3)
        % (c1'*x + )^alpha*(c2'x+ )^(1-alpha) - norm(c3*x + )
        % x3 - x2exp(x1/x2) >= 0
        % d/dx = [-exp(x1) -(exp(x1/x2)-x1/x2*exp(x1/x2)) 1]  
        if x(1)*x(2)==0
            a = 1;b = 1;
        else
            a = x1^(alpha-1)*x2^(1-alpha);
            b = x1^alpha*x2^(-alpha);
        end
        if nargin == 2
            c1=c1(model.linearindicies);
            c2=c2(model.linearindicies);
            c3=c3(:,model.linearindicies);
            smoothed = sqrt(10^-15 + x3'*x3);            
            temp = c1*alpha*a + (1-alpha)*c2*b - (x3'*c3)/smoothed;
            conederiv = [conederiv -temp'];
        else
            smoothed = sqrt(10^-15 + x3'*x3);
            temp = c1*alpha*a + (1-alpha)*c2*b - (x3'*c3)/smoothed;                       
            conederiv = [conederiv -(temp*dzdx)'];           
        end
        top = top + n;
    end
end
conederiv = conederiv';

function [dg,reordering] = computeSDPDeriv(model,xevaled,newdxx)
top = startofSDPCone(model.K);
newF = [];newcuts = 1;
B=model.F_struc(:,2:end)*newdxx;
global sdpLayer
reordering = [];

for i = 1:length(model.K.s)
    n = model.K.s(i);
    X = model.F_struc(top:top+n^2-1,:)*[1;xevaled];
    X = full(reshape(X,n,n));
    [d,v] = eig(X);
    if ~isempty(sdpLayer.eigenVectors)
        use = zeros(1,n);
        used = zeros(1,size(sdpLayer.eigenVectors{i},2));
        for j = 1:n
            if ~isempty(sdpLayer.eigenVectors{i})
                Xv = X*sdpLayer.eigenVectors{i};
            end
            for k = 1:size(sdpLayer.eigenVectors{i},2)
                s = norm(Xv(:,k)-v(j,j)*sdpLayer.eigenVectors{i}(:,k));
                if s <= 1e-12
                    if ~used(k) && ~use(j)
                        use(j) = k;
                        used(k) = 1;
                    end
                end
            end
        end
        sdpLayer.eigenVectors{i} = [sdpLayer.eigenVectors{i} d];
        if size(sdpLayer.eigenVectors{i},2) > 10*size(sdpLayer.eigenVectors{i},1)
            sdpLayer.eigenVectors{i} = sdpLayer.eigenVectors{i}(:,end-10*size(sdpLayer.eigenVectors{i},1):end);
        end        
    else
        sdpLayer.eigenVectors{i} = d;
    end
    %  use
%    if any( find(use))
%        disp(use)
%    end
    d(:,find(use)) = sdpLayer.eigenVectors{i}(:,use(find(use)));
    %   if size(sdpLayer.eigenVectors{i},2) > 5*n
    %       sdpLayer.eigenVectors{i} = sdpLayer.eigenVectors{i}(:,end-5*n:end);
    %   end
    newSDPblock = [];
    nZeroVectors = length(find(abs(diag(v))<=1e-8));
    try
        if any(nZeroVectors)
            if ~isempty(sdpLayer.nullVectors{i})
                tt = [];
                for ee = 1:size(sdpLayer.nullVectors{i},2)
                    tt = [tt sdpLayer.nullVectors{i}(:,ee)'*X*sdpLayer.nullVectors{i}(:,ee)];
                end
                notok = find(tt >= 2e-6);
                if any(notok)
                   % 'Prune'
                    sdpLayer.nullVectors{i}(:,find(notok))=[];
                end
            end
            if isempty(sdpLayer.nullVectors{i})
                sdpLayer.nullVectors{i} = d(:,1:nZeroVectors);
                %'New batch'
            elseif  nZeroVectors > size(sdpLayer.nullVectors{i},2)
                %'Expanding'
                S = [sdpLayer.nullVectors{i} d(:,nZeroVectors+1:end)];
                S = null(S');
                %S = null(sdpLayer.nullVectors{i}');
                sdpLayer.nullVectors{i} = [sdpLayer.nullVectors{i} S(:,1:(nZeroVectors-size(sdpLayer.nullVectors{i},2)))];
            elseif nZeroVectors < size(sdpLayer.nullVectors{i},2)
                %'Contracting'
            end
            %  'Using'
            d(:,1:nZeroVectors) = sdpLayer.nullVectors{i}(:,1:nZeroVectors);
        end
    catch
        %  'Fail'
    end
    pre_reordering = 1:n;
    %     if any(abs(diff(diag(v))) <=1e-6)
    %         groupStart = 1;
    %         vv = diag(v);
    %         while groupStart < n
    %             thisGroup = find(abs(vv(groupStart+1:end)-vv(groupStart))<=1e-6);
    %             if any(thisGroup)
    %                 members = [groupStart groupStart+thisGroup(:)'];
    %                 [~,loc] = sort((1:n)*abs(d(:,members)));
    %                 pre_reordering(members) = members(loc);
    %                 groupStart = groupStart + length(members);
    %             else
    %                 groupStart = groupStart + 1;
    %             end
    %         end
    %     end
    for m = 1:min(sdpLayer.n,model.K.s(i))
        newrow = [];
        newrow = oneSDPGradient(d(:,m),B(top:top+n^2-1,:));
        newSDPblock = [newSDPblock;newrow];
        newcuts = newcuts + 1;
    end
    if strcmpi(model.options.slayer.algorithm,'convex') || isequal(model.options.slayer.algorithm,1)
        newSDPblock = cumsum(newSDPblock);
    end
    if ~isempty(sdpLayer.oldGradient{i})
        [reordering] = matchGradientRows(newSDPblock,sdpLayer.oldGradient{i});
        newSDPblock = newSDPblock(reordering,:);
        reordering = reordering(pre_reordering);
    end
    %  newSDPblock
    sdpLayer.oldGradient{i} = newSDPblock;
    sdpLayer.reordering{i} = reordering;
    newF = [newF;newSDPblock];
    top = top + n^2;
end
dg = -newF;

function newdxx = computeMonomialVariableJacobian(model,x)
n = length(model.c);
linearindicies = model.linearindicies;

news = model.fastdiff.news;
allDerivemt = model.fastdiff.allDerivemt;
c = model.fastdiff.c;

if model.fastdiff.univariateDifferentiates
    zzz = c.*(x(model.fastdiff.univariateDiffMonom).^model.fastdiff.univariateDiffPower);
else
    %  X = repmat(x(:)',length(c),1);
    O = ones(length(c),length(x));
    nz = find(allDerivemt);
    %  O(nz) = X(nz).^allDerivemt(nz);
    O(nz) = x(ceil(nz/length(c))).^allDerivemt(nz);
    zzz = c.*prod(O,2);
end
newdxx = model.fastdiff.newdxx;
newdxx(model.fastdiff.linear_in_newdxx) = zzz;

function r = matchGradientRows(X,Y)
C = distanceMatrix(X,Y);
[r,m,u] = matchpairs(C,1e12);
r = r(:,1);

function C = distanceMatrix(X,Y)
C = zeros(size(X,1));
for i = 1:size(X,1)
    for j = 1:size(X,1)
        C(i,j) = norm(X(i,:)-Y(j,:),1);               
    end
end

function  g = reorderEigenvalueG(g,model)
global sdpLayer
if ~isempty(sdpLayer.reordering{1})
    g_nonsdp = g(1:end-sum(min(sdpLayer.n,model.K.s)));
    g_sdp = g(end-sum(min(sdpLayer.n,model.K.s))+1:end);
    reordering = [];
    top = 0;
    for i = 1:length(model.K.s)
        reordering = [reordering;top+sdpLayer.reordering{i}];
        top = top + min(sdpLayer.n,model.K.s(i));
    end
    g_sdp = g_sdp(reordering);
    g = [g_nonsdp;g_sdp];
end

function newrow = oneSDPGradient(d,B,top)
n = size(d,1);
newrow = [];
for j = 1:size(B,2)
    newrow = [newrow d'*reshape(B(:,j),n,n)*d];
end
