function output = callpowersolver(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
K       = model.K;
ub      = model.ub;
lb      = model.lb;

% Create the parameter structure
pars = options.powersolver;

% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

if options.savedebug
    save powersolverdebug model
end

% *********************************************
% Call SeDuMi
% *********************************************
if options.showprogress;
    showprogress(['Calling ' model.solver.tag],options.showprogress);
end

K.p = K.p';
K.e = 0;
K.pd = 'd';
solvertime = tic;
[x_s,y_s,info] = powersolver(-F_struc(:,2:end)',-c,F_struc(:,1),K);
solvertime = toc(solvertime);

% Internal format
Primal = y_s; 
Dual   = x_s;

problem = 0;

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