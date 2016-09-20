function output = callsdpt331(interfacedata)

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

% Convert from internal (sedumi-like) format
[blk,A,C,b,oldKs]=sedumi2sdpt3(F_struc(:,1),F_struc(:,2:end),c,K,options.sdpt3.smallblkdim);

options.sdpt3.printyes=double(options.verbose);
options.sdpt3.expon=options.sdpt3.expon(1);
if options.savedebug
    ops = options.sdpt3;
    save sdpt3debug blk A C b ops x0
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = tic;
if options.verbose==0 % SDPT3 does not run silent despite printyes=0!
   evalc('[obj,X,y,Z,info,runhist] =  sqlp(blk,A,C,b,options.sdpt3,[],x0,[]);');
else
    [obj,X,y,Z,info,runhist] =  sqlp(blk,A,C,b,options.sdpt3,[],x0,[]);    
end
solvertime = toc(solvertime);

% Create YALMIP dual variable and slack
Dual = [];
Slack = [];
top = 1;
if K.f>0
    Dual = [Dual;X{top}(:)];
    Slack = [Slack;Z{top}(:)];
    top = top+1;
end
if K.l>0
    Dual = [Dual;X{top}(:)];
    Slack = [Slack;Z{top}(:)];
    top = top + 1;
end
if K.q(1)>0
    Dual = [Dual;X{top}(:)];
    Slack = [Slack;Z{top}(:)];
    top = top + 1;
end
if K.s(1)>0  
    % Messy format in SDPT3 to block and sort small SDPs
    u = blk(:,1);
    u = find([u{:}]=='s');
    s = 1;
    for top = u
        ns = blk(top,2);ns = ns{1};
        k = 1;
        for i = 1:length(ns)
            Xi{oldKs(s)} = X{top}(k:k+ns(i)-1,k:k+ns(i)-1);
            Zi{oldKs(s)} = Z{top}(k:k+ns(i)-1,k:k+ns(i)-1);
            s = s + 1;                 
            k = k+ns(i);
        end
    end 
    for i = 1:length(Xi)
        Dual = [Dual;Xi{i}(:)];     
        Slack = [Slack;Zi{i}(:)];     
    end
end
Primal = -y;  % Primal variable in YALMIP

% Convert error code
switch info.termcode
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