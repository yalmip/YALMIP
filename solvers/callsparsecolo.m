function output = callsparsecolo(interfacedata)

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

A = -F_struc(:,2:end)';
C = F_struc(:,1);
b = -c;

if options.savedebug
    ops = options.sparsecolo;
    save sparsecolodebug K A C b
end

ops = options.sparsecolo;
ops.SDPsolver = lower(interfacedata.solver.sdpsolver.tag);
ops.SDPAoptions = options.sdpa;
ops.sedumipar = options.sedumi;
ops.sdpt3OPTIONS = options.sdpt3;
ops.sdpt3OPTIONS.printlevel = options.verbose;
ops.printlevel = options.verbose;
if options.verbose==0
    ops.SDPAoptions.print = 'no';
else
    ops.SDPAoptions.print = 'display';
end

if options.savedebug
    save sparsecolodebug K A C b ops
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = tic;
if options.verbose==0 % Sparsecolo does not run silent
    evalc('[x,y,infoCoLO,cliqueDomain,cliqueRange,LOP] = sparseCoLO(A,b,C,K,[],ops);');
else    
    [x,y,infoCoLO,cliqueDomain,cliqueRange,LOP] = sparseCoLO(A,b,C,K,[],ops);
end
solvertime = toc(solvertime);

% Create YALMIP dual variable and slack
Dual = x;
Primal = y;


switch  lower(interfacedata.solver.sdpsolver.tag)
    case 'sdpt3'
        switch infoCoLO.SDPsolver.termcode
            case 0
                problem = 0; % No problems detected
            case {-1,-5}
                problem = 5; % Lack of progress
            case {-2,-3,-4,-7}
                problem = 4; % Numerical problems
            case -6
                problem = 3; % Maximum iterations exceeded
            case -10
                problem = 7; % YALMIP sent incorrect input to solver
            case 1
                problem = 2; % Dual feasibility
            case 2
                problem = 1; % Primal infeasibility
            otherwise
                problem = -1; % Unknown error
        end
    case 'sedumi'
        problem = sedumicode(infoCoLO.SDPsolver,options);
        
    case 'sdpa'
        switch (infoCoLO.SDPsolver.phasevalue)
            case 'pdOPT'
                problem = 0;
            case {'noINFO','pFEAS','dFEAS','pdFEAS'}
                problem = 3;
            case 'pFEAS_dINF'
                problem = 2;
            case 'pINF_dFEAS'
                problem = 1;
            case 'pUNBD'
                problem = 2;
            case 'dUNBD'
                problem = 1;
            case 'pdINF'
                problem = 12;
            otherwise
                problem = -1;
        end
        
end
infostr = yalmiperror(problem,interfacedata.solver.tag);

if options.savesolveroutput
    solveroutput.obj = obj;
    solveroutput.X = X;
    solveroutput.y = y;
    solveroutput.Z = Z;
    solveroutput.info = info;
    solveroutput.runhist = runhist;
else
    solveroutput = [];
end

if options.savesolverinput
    solverinput.blk = blk;
    solverinput.A   = A;
    solverinput.C   = C;
    solverinput.b   = b;
    solverinput.X0   = [];
    solverinput.y0   = x0;
    solverinput.Z0   = [];
    solverinput.options   = options.sdpt3;
else
    solverinput = [];
end

% Standard interface 
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);


function problem = sedumicode(info,options)
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

