function output = calloptidsdp(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

model = yalmip2opticsdp(interfacedata);
if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end

if options.savedebug
    save csdpdebug model
end

solvertime = tic;
model.ops = csdpset;
[y,fvals,exitflag,stats,X] = csdp(model.f,model.A,model.b,model.lb,model.ub,model.sdcone,model.y0,model.ops);
solvertime = toc(solvertime);

% Create dual variable in internal format
if options.saveduals
    if K.l > 0
        top = 1;
        D_struc = X{1};
    else
        top = 0;
        D_struc = [];
    end
    if K.s(1) > 0
        for j = 1:length(K.s)                      
            D_struc = [D_struc;X{j+top}(:)];
        end
    end
else
    D_struc = [];
end

x = y;  % Our notation do not coincide ...
switch exitflag
    case 0
       problem = 0;
    case 1
        problem = 2;
    case 2
        problem = 1;
    case {3,5,6,7,8,9,10}
        problem = 4;
    case {4,-27}
        problem = 3;
    case -50
        problem = 16;
    case 11
        problem = 7;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,interfacedata.solver.tag);

if options.savesolveroutput   
	solveroutput.y = y;
    solveroutput.fvals = fvals;
    solveroutput.exitflag = exitflag;
    solveroutput.stats = stats;
    solveroutput.X = X;
else
	solveroutput = [];
end

if options.savesolverinput
	solverinput.model = model;	
else
	solverinput = [];
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);