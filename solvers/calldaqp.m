function output = calldaqp(interfacedata)

% Get problem data
options = interfacedata.options;
H       = full(2*interfacedata.Q);
f       = full(interfacedata.c);
F_struc = interfacedata.F_struc;
K       = interfacedata.K;
lb      = interfacedata.lb;
ub      = interfacedata.ub;
bin_vars= interfacedata.binary_variables;
% Extract constraints
if ~isempty(F_struc)
  A = full(-F_struc(:,2:end));
  b = full(F_struc(:,1));
  meq = K.f;
  m = length(b)-meq;
else
    A = []; b = [];
	m=0; meq=0;
end
bupper = full([ub;b]);
blower = full([lb;b(1:meq);-inf(m,1)]);
% Set constraint type
sense =  [zeros(length(ub),1,'int32');5*ones(meq,1,'int32');zeros(m,1,'int32')];
sense(find(isinf(lb)&isinf(ub))) = int32(4); % "ignore" if lb and ub are inf
sense(bin_vars) = int32(16);

if options.savedebug
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end

% Solve with DAQP 
d = daqp();
d.settings(options.daqp);
solvertime = tic;
exitflag = d.setup(H,f,A,bupper,blower,sense);
if(exitflag > 0)
  [xstar,fval,exitflag,info] = d.solve();
else
  xstar = []; info= [];
end
solvertime = toc(solvertime);

switch exitflag 
    case 1
        problem = 0;
    case 2
        problem = 0;
    case -1
        problem = 1;
    case -2
        problem = 5;
    case -3 
        problem = 2;
    case -4
        problem = 3;
    case -5 
        problem = 14;
    case -6
        problem = 6;
    otherwise
        problem = 9;
end

% Standard interface
Primal      = xstar;
Dual        = info.lambda(length(ub)+1:end);
infostr     = yalmiperror(problem,interfacedata.solver.tag);
if ~options.savesolverinput
    solverinput = [];
else
    solverinput = model;
end
if ~options.savesolveroutput
    solveroutput = [];
else
    solveroutput = info;
end

% Standard interface
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);
