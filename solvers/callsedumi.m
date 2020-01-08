function output = callsedumi(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
K       = model.K;
ub      = model.ub;
lb      = model.lb;

% Create the parameter structure
pars = options.sedumi;
pars.fid = double(options.verbose);

% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

% Avoid bug (by design?) in SeDuMi on 1D second order cones
if any(K.q == 2)
    [F_struc,K] = fix1Dqcone(F_struc,K);
end

if options.savedebug
    save sedumidebug F_struc c K pars
end

% *********************************************
% Call SeDuMi
% *********************************************
if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end

problem = 0;  
warnState = warning;
solvertime = tic;
try
    [x_s,y_s,info] = sedumi(-F_struc(:,2:end),-c,F_struc(:,1),K,pars);
catch
    if strfind(lasterr,'Out of memory')
        error(lasterr)
    end
    try
        if options.verbose > 0
            disp(' ');
            disp('SeDuMi had unexplained problems, maybe due to linear dependence?')
            disp('YALMIP tweaks the problem (adds 1e6 magnitude bounds on all variables) and restarts...')            
            disp(' ');
        end
        % Boring issue in sedumi for trivial problem min x+y, s.t x+y>0
        n = length(c);
        [F_struc,K] = addStructureBounds(F_struc,K,ones(n,1)*1e6,-ones(n,1)*1e6);
        [x_s,y_s,info] = sedumi(-F_struc(:,2:end),-c,F_struc(:,1),K,pars);
        x_s((1:2*n)+K.f)=[];
        K.l=K.l-2*n;
    catch
        disp('Nope, unexplained crash in SeDuMi! (could be memory issues or wrong binary)')
        disp('Make sure you have a recent and compiled version')
        disp(' ')        
        disp(' ')        
        disp('For better diagnostics, use sdpsettings(''debug'',1)')        
        disp(' ')        
        disp(' ')                
    end
end
warning(warnState);
solvertime = toc(solvertime);



% Internal format
Primal = y_s; 
Dual   = x_s;

temp = info.pinf;
pinf = info.dinf;
dinf = temp;

% Check for reported errors
if (pinf==0) &  (dinf==0)
    problem = 0; % No problems
end

% We can only report one error, use priorities
if (problem==0) & (pinf==1)
    problem = 1; % Primal infeasability
end

if (problem==0) & (dinf==1)
    problem = 2; % Dual infeasability
end

if (problem==0) & (info.numerr==1) | (info.numerr==2)
    problem = 4; %Numerical problems
end

if (problem==0) & (info.iter >= options.sedumi.maxiter)
    % Did we need exactly maxiter iterations to find optimum
    if (pinf==0) & (dinf==0) & (info.numerr==0)
        problem = 0; % Yes
    else
        problem = 3; % No, we are not optimal yet
    end
end

if (problem==0) & (info.feasratio<0.98) 
    problem = 4; 
end

% Fix for cases not really explained in documentation of sedumi?
if (abs(info.feasratio+1)<0.1) & (pinf==0) &  (dinf==0)
    problem = 1; 
end
if (abs(info.feasratio+1)<0.1) & (pinf==0) &  (dinf==0) & (c'*y_s<-1e10)
    problem = 2; 
end


infostr = yalmiperror(problem,model.solver.tag);

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.A = -F_struc(:,2:end);
    solverinput.c = F_struc(:,1);
    solverinput.b = -c;
    solverinput.K = K;
    solverinput.pars = pars;
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
    solveroutput.x = x_s;
    solveroutput.y = y_s;
    solveroutput.info = info;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);

 
function [new_F_struc,K] = fix1Dqcone(F_struc,K); 

new_F_struc = F_struc(1:(K.l+K.f),:);
top = K.f + K.l + 1;
for i = 1:length(K.q)
    if K.q(i) == 2
      new_F_struc = [new_F_struc;F_struc(top,:);F_struc(top+1,:);F_struc(top,:)*0];
      K.q(i) = 3;
      top = top + 2;
    else
      new_F_struc = [new_F_struc;F_struc(top:top+K.q(i)-1,:)];
      top = top + K.q(i);
    end
end
new_F_struc = [new_F_struc;F_struc(top:end,:)];
