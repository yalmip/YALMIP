function output = calldsdp(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

% Bounded variables converted to constraints
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

% Convert from SeDuMi format
[b,AC] = sedumi2dsdp5(F_struc,c,K);

options.dsdp.printyes = (options.verbose>0);
options.dsdp.print = (options.verbose>0);

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end

if options.savedebug
    pars = options.dsdp;
    save dsdpdebug AC b pars
end

solvertime = tic;
if isempty(x0)
    if options.saveduals | options.dimacs | options.savesolveroutput
        if options.verbose==0 % to fix display bug reported from user
            evalc('[STAT,y,X] =  dsdp(b,AC,options.dsdp);');
        else
            [STAT,y,X] =  dsdp(b,AC,options.dsdp);
        end
    else
        if options.verbose==0 % to fix display bug reported from user
            evalc('[STAT,y] =  dsdp(b,AC,options.dsdp);');
        else
            [STAT,y] =  dsdp(b,AC,options.dsdp);;
        end
    end
else
    if options.saveduals | options.dimacs | options.savesolveroutput
        if options.verbose==0 % to fix display bug reported from user
            evalc('[STAT,y,X] =  dsdp(b,AC,options.dsdp,x0);');
        else
            [STAT,y,X] =  dsdp(b,AC,options.dsdp,x0);
        end
    else
        if options.verbose==0 % to fix display bug reported from user
            evalc('[STAT,y] =  dsdp(b,AC,options.dsdp,x0);');
        else
            [STAT,y] =  dsdp(b,AC,options.dsdp,x0);
        end
    end
end
solvertime = toc(solvertime);

% Create dual variable in internal format
if options.saveduals | options.dimacs
    Dual = [];
    for i = 1:length(X)
        if isequal('SDP',AC{i,1})
            XX = dmat(X{i});
        else
            XX = X{i};
        end
        Dual = [Dual;XX(:)];
    end
else
    Dual = [];
end

Primal = y;  % Our notation do not coincide ...

% Backwards compatibility
if isfield(STAT,'dual')
    switch STAT.dual
        case 'Infeasible'
            STAT.termcode = 2;
        case 'Unbounded'
            STAT.termcode = 1;
        otherwise
    end
end

switch STAT.termcode
    case 0
        if STAT.iterates>=options.dsdp.maxit
            problem = 3;
        else
            problem = 0;
        end
    case 1
        problem = 2;
    case 2
        problem = 1;
    case -3
        problem = 4;
    otherwise
        problem = -1;
end

if options.savesolveroutput
	solveroutput.STAT = STAT;
	solveroutput.X = X;
	solveroutput.y = y;
else
	solveroutput = [];
end

if options.savesolverinput
	solverinput.AC = AC;	
	solverinput.b = b;
	solverinput.pars = options.dsdp;
else
	solverinput = [];
end

% Standard interface 
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);