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
[C,A,b,blk] = sedumi2dsdp(F_struc,c,K);

% Quadratic not yet implemented?
options.dsdp.dual_quadratic=spalloc(length(c),length(c),0);

options.dsdp.printyes = (options.verbose>0);
if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end

if options.savedebug
    save dsdpdebug A C b options.dsdp
end

solvertime = tic;
if isempty(x0)
    if options.saveduals | options.savesolveroutput
        if options.verbose==0 % to fix display bug reported from user
            evalc('[STAT,y,X] =  dsdp(A,C,b,options.dsdp);');
        else
            [STAT,y,X] =  dsdp(A,C,b,options.dsdp);
        end
    else
        if options.verbose==0 % to fix display bug reported from user
            evalc('[STAT,y] =  dsdp(A,C,b,options.dsdp);');
        else
            evalc('[STAT,y] =  dsdp(A,C,b,options.dsdp);');
        end
    end
else
    if options.saveduals | options.savesolveroutput
        if options.verbose==0 % to fix display bug reported from user
            evalc('[STAT,y,X] =  dsdp(A,C,b,options.dsdp,x0);');
        else
            [STAT,y,X] =  dsdp(A,C,b,options.dsdp,x0);
        end
    else
        if options.verbose==0 % to fix display bug reported from user
            evalc('[STAT,y] =  dsdp(A,C,b,options.dsdp,x0);');
        else
            [STAT,y] =  dsdp(A,C,b,options.dsdp,x0);
        end
    end
end
solvertime = toc(solvertime);

% Create dual variable in internal format
if options.saveduals
    D_struc = [];
    for i = 1:length(X)
        D_struc = [D_struc;X{i}(:)];
    end
else
    D_struc = [];
end

x = y;  % Our notation do not coincide ...
switch STAT.termcode
case 0 
	if STAT.iterates>=options.dsdp.maxit
		problem = 3;
	else
		problem = 0;
	end
%     if STAT.gaphist (end) > options.dsdp.gaptol
%         problem = 5;
%     end
case 1
	problem = 2;
case 2
	problem = 1;
case -3
    problem = 3;
otherwise
	problem = -1;
end
infostr = yalmiperror(problem,interfacedata.solver.tag);

if options.savesolveroutput
	solveroutput.STAT = STAT;
	solveroutput.X = X;
	solveroutput.y = y;
else
	solveroutput = [];
end

if options.savesolverinput
	solverinput.A = A;
	solverinput.C = C;
	solverinput.b = c;
	solverinput.pars = options.dsdp;
else
	solverinput = [];
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);