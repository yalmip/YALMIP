function output = callvsdp(interfacedata)

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

% Convert from internal (sedumi-like) format to VSDPs sdpt3-like format
[blk,A,C,b]=sedumi2vsdp(F_struc(:,1),F_struc(:,2:end),c,K);

if options.savedebug
    ops = options.sdpt3;
    save sdpt3debug blk A C b ops -v6
end

% Solver to be used in VSDP
options.vsdp.model = interfacedata;

solvertime = tic;     
[objt,Xt,yt,Zt,info] = mysdps_yalmip(blk,A',C,b,options);

% Compute rigorous lower bound (default behaviour)
if options.vsdp.verifiedlower
    [fL, Y, dL] = vsdplow_yalmip(blk,A',C,b,Xt,yt,Zt,[],options);
    if isnan(Y)
        info(1) = 11;
    end
    %[fL, Y, dL] = vsdplow(blk,A,C,b,Xt,yt,Zt)    
else
    Y = [];
    fL = [];
    dL = [];
end

% Compute rigorous lower bound
if options.vsdp.verifiedupper
    %[fL, Y, dL] = vsdplow_yalmip(blk,A,C,b,[],[],[],[],options)    
    % Now compute rigorous lower bound
    [fU, X, lb] = vsdpup_yalmip(blk,A',C,b,Xt,yt,Zt,[],options);
    %[fU, X, lb] = vsdpup(blk,A,C,b,Xt,yt,Zt);    
else
    fU = [];
    X = [];
    lb = [];
end

solvertime = toc(solvertime);

Dual = [];
Slack = [];
top = 1; 
if K.f>0
    Dual = [Dual;Xt{top}(:)];
    Slack = [Slack;Zt{top}(:)];
    top = top+1;
end
if K.l>0
    Dual = [Dual;[Xt{1:K.l}]'];
    Slack = [Slack;[Zt{1:K.l}]'];
    top = top + K.l;
end
if K.q(1)>0
    Dual = [Dual;Xt{top}(:)];
    Slack = [Slack;Zt{top}(:)];
    top = top + 1;
end
if K.s(1)>0     
    for i = 1:length(K.s)
        Dual = [Dual;Xt{top+i-1}(:)];     
        Slack = [Slack;Zt{top+i-1}(:)];     
    end
end

if ~isempty(Y) & ~isnan(Y)
    Primal = Y;  % Primal variable in YALMIP
else
    Primal = yt;  % Primal variable in YALMIP
end

% Convert error code
switch info(1)
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
    case 11
        problem = 11;
    otherwise
        problem = -1; % Unknown error
end
infostr = yalmiperror(problem,interfacedata.solver.tag);

% always save output
if options.savesolveroutput
    solveroutput.objt = objt;
    solveroutput.Xt = Xt;
    solveroutput.yt = yt;
    solveroutput.Zt = Zt;
    solveroutput.fL = fL;
    solveroutput.Y  = Y;
    solveroutput.dL  = dL;    
    solveroutput.fU = fU;
    solveroutput.X  = X;
    solveroutput.lb  = lb;
    solveroutput.info = info;
 else
    solveroutput = [];
end

if options.savesolverinput
    solverinput.blk = blk;
    solverinput.A   = A;
    solverinput.C   = C;
    solverinput.b   = b; 
    solverinput.options   = options.sdpt3;
else
    solverinput = [];
end

% Standard interface 
output = createOutputStructure(Primal,Dual,Slack,problem,infostr,solverinput,solveroutput,solvertime);


function [blk,A,C,b]=sedumi2vsdp(Cin,Ain,c,K);

C = {};
A = {};
b  = -c;
blk = {};
top = 1;
k = 1;

if K.f > 0
    C{k,1} = Cin(top:top+K.f-1);
%    A{k} = -Ain(top:top+K.f-1,:);
    for j = 1:length(c)
        A{j,k} = -Ain(top:top+K.f-1,j);
    end
    blk{k,1} = 'u';
    blk{k,2} = K.f;
    top = top + K.f;
    k = k + 1;
end

if K.l > 0
    K.s = [repmat(1,1,K.l) K.s];
    K.s(K.s == 0) = [];
%     C{k,1} = Cin(top:top+K.l-1);
% %    A{k} = -Ain(top:top+K.l-1,:);
%     for j = 1:length(c)
%         A{j,k} = -Ain(top:top+K.l-1,j);
%     end
%     blk{k,1} = 'l';
%     blk{k,2} = K.l;
%     top = top + K.l;
%     k = k + 1;
end

if K.s(1)>0
    for i = 1:length(K.s)
        C{k,1} = reshape(Cin(top:top+K.s(i)^2-1),K.s(i),K.s(i));
        for j = 1:length(c)
            A{j,k} = reshape(-Ain(top:top+K.s(i)^2-1,j),K.s(i),K.s(i));
        end
        blk{k,1} = 's';
        blk{k,2} = K.s(i);
        k = k + 1;
        top = top + K.s(i)^2;
    end
end

A = A';