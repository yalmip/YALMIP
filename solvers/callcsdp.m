function output = callcsdp(interfacedata);

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
pars   = interfacedata.options.csdp;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

% Bounded variables converted to constraints
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

pars.printlevel=options.verbose;

if options.savedebug
    save csdpdebug F_struc c K pars
end
    
if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = tic;
if options.verbose==0 % to fix display bug reported from user
    evalc('[x_s,y_s,z_s,info]=csdp(-F_struc(:,2:end),-c,F_struc(:,1),K,pars);');
else
    [x_s,y_s,z_s,info]=csdp(-F_struc(:,2:end),-full(c),F_struc(:,1),K,pars);
end
solvertime = toc(solvertime);

% We solve dual problem with CSDP
Dual = x_s;
x    = y_s;

switch info
case 0
    problem = 0;
case 1
    problem = 2;
case 2
    problem = 1;
case {3,5,6,7,8,9}
    problem = 4;
case 4
    problem = 3;
otherwise
    problem = 9;
end 
infostr = yalmiperror(problem,interfacedata.solver.tag);

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
output = createOutputStructure(x,Dual,z_s,problem,infostr,solverinput,solveroutput,solvertime);