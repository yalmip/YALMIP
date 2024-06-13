function output = callcdcs(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
K       = model.K;
ub      = model.ub;
lb      = model.lb;


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

options.cdcs.verbose = options.verbose;
aK.f = K.f;
aK.l = K.l;
aK.q = K.q;
aK.s = K.s;
K = aK;

% Create the options structure for cdcs
opts = options.cdcs;
opts.verbose = double(options.verbose);

if options.savedebug
    save cdcsdebug F_struc c K opts
end

% *********************************************
% Call CDCS
% *********************************************
if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end
 
% GF: this shoud now work. Outputs x_s, y_s are in standard SeDuMi format
problem = 0;  
solvertime = tic;
try
    [x_s,y_s,z_s,info] = cdcs(-F_struc(:,2:end),-c,F_struc(:,1),K,options.cdcs);
    % Internal format
    Primal = y_s; 
    Dual   = x_s;
catch
    disp('Unexpected crash in CDCS!')
    disp('Make sure you have a recent working version of CDCS (test with CDCSTest)')
    problem = 9;
end
solvertime = toc(solvertime);

% Set problem code from solver if successful run
if problem~=9
    if info.problem<4
        % CDCS problem code matches YALMIP
        % 0: successfully solved
        % 1: primal infeasible
        % 2: dual infeasible (unbounded objective)
        % 3: max number of iterations reached
        problem = info.problem;
    else
        % Solution probably good, but error in post-processing routine
        % Please check if problem code is appropriate
        problem = 11;
    end
end

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.Atveriny = -F_struc(:,2:end);
    solverinput.c = F_struc(:,1);
    solverinput.b = -c;
    solverinput.K = K;
    solverinput.options = options.cdcs;
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
    solveroutput.x = x_s;
    solveroutput.y = y_s;
    solveroutput.z = z_s;
    solveroutput.info = info;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(Primal,Dual,[],problem,model.solver.tag,solverinput,solveroutput,solvertime);

 
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