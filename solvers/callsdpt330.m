function output = callsdpt330(interfacedata)

% Author Johan Löfberg
% $Id: callsdpt330.m,v 1.3 2005-05-07 13:53:20 joloef Exp $ 

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
    [F_struc,K] = addbounds(F_struc,K,ub,lb);
end

% Convert from internal (sedumi) format
[C,A,b,blk] = sdpt3struct2sdpt3block(F_struc,c,K);

% To be sure
global CACHE_SIZE   % cache size in kbytes 
global LOOP_LEVEL   % loop unrolling level
CACHE_SIZE = options.sdpt3.CACHE_SIZE; 
LOOP_LEVEL = options.sdpt3.LOOP_LEVEL;

solvertime = clock;
A = svec(blk,A,ones(size(blk,1),1));

if options.savedebug
    save sdpt3debug blk A C b x0 options.sdpt3
end

showprogress('Calling SDPT3',options.showprogress);
if options.verbose==0
    evalc('[obj,X,y,Z,gaphist,infeashist,info,Xiter,yiter,Ziter] =  sqlp(blk,A,C,b,[],x0,[],options.sdpt3);');
else
    [obj,X,y,Z,gaphist,infeashist,info,Xiter,yiter,Ziter] =  sqlp(blk,A,C,b,[],x0,[],options.sdpt3);
end
solvertime = etime(clock,solvertime);

% Create dual variable in internal format
D_struc = [];
if K.l>0
    D_struc = [D_struc;X{1}(:)];
    if length(X)>1
        X = {X{2:end}};
    end
end
if K.q(1)>0
    D_struc = [D_struc;X{1}(:)];
    if length(X)>1
        X = {X{2:end}};
    end
end
if K.r(1)>0
    D_struc = [D_struc;X{1}(:)];
    if length(X)>1
        X = {X{2:end}};
    end
end
if K.s(1)>0
    top = 1;
    j = 1;
    for i = 1:1:length(K.s)
        Xi = X{j}(top:top+K.s(i)-1,top:top+K.s(i)-1);
        D_struc = [D_struc;Xi(:)];
        top = top + K.s(i);
        if top>length(X{j})
            j = j + 1;
            top = 1;
        end
    end
end

x = y;  % Our notation do not coincide ...

% Convert error code
switch info(1)
    case 0
        problem = 0; % No problems detected
    case -1 
        problem = 5; % Lack of progress
    case {-2,-3,-4,-5}
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
    solveroutput.gaphist = gaphist;
    solveroutput.infeashist = infeashist;
    solveroutput.info = info;
    solveroutput.Xiter = Xiter;
    solveroutput.yiter = yiter;
    solveroutput.Ziter = Ziter;
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
output.Primal      = x(:);
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;