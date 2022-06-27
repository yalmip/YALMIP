function output = calldaqp(interfacedata)

options = interfacedata.options;
model = yalmip2quadprog(interfacedata);

if options.savedebug
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end


% Define QP
[m,n]= size(model.A);
meq = length(model.beq);
H = full(model.Q);
f = full(model.c);
A = full([model.A;model.Aeq]);
if(all(isinf(model.lb))&all(isinf(model.ub)))
  lb=[];ub=[];
else 
  lb = model.lb; ub=model.ub;
end
bupper = full([ub;model.b;model.beq]);
blower = full([lb;-inf(m,1);model.beq]);
sense =  [zeros(length(ub)+m,1,'int32');5*ones(meq,1,'int32')];

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
Dual        = []; % TODO: add dual information 
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
