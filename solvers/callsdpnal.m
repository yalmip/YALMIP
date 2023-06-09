function output = callsdpnal(interfacedata)

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

[blk,A,C,b,oldKs]=sedumi2sdpt3(F_struc(:,1),F_struc(:,2:end),c,K,1);

if options.savedebug
    ops = options.sdpnal;
    save sdpnaldebug blk A C b ops -v6
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = tic;
if options.verbose==0
    evalc('[obj,X,s,y,Z,Z2,y2,v,info,runhist] = sdpnalplus(blk,A,C,b,[],[],[],[],[],options.sdpnal);');
else
    [obj,X,s,y,Z,Z2,y2,v,info,runhist] = sdpnalplus(blk,A,C,b,[],[],[],[],[],options.sdpnal);          
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
if any(K.q)
    Dual = [Dual;X{top}(:)];
    Slack = [Slack;Z{top}(:)];
    top = top + 1;
end
if any(K.s)  
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

% No error code available
if isfield(info,'termcode')
    switch info.termcode
        case -1
            problem = 4;
        case -2
            problem = 3;
        case 0
            problem = 0;
        case 1
            problem = 5;
        case 2
            problem = 3;
        otherwise
            problem = 11;
    end
else
    if isfield(info,'msg')
        if isequal(info.msg,'maximum iteration reached')
            problem = 3; 
        else
            problem = 9;
        end
    end
end

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
    solverinput.options   = options.sdpnal;
else
    solverinput = [];
end

% Standard interface 
output = createOutputStructure(Primal,Dual,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);


function [F_struc,K] = deblock(F_struc,K);
X = any(F_struc(end-K.s(end)^2+1:end,:),2);
X = reshape(X,K.s(end),K.s(end));
[v,dummy,r,dummy2]=dmperm(X);
blks = diff(r);

lint = F_struc(1:end-K.s(end)^2,:);
logt = F_struc(end-K.s(end)^2+1:end,:);

newlogt = [];
for i = 1:size(logt,2)
    temp = reshape(logt(:,i),K.s(end),K.s(end));
    temp = temp(v,v);
    newlogt = [newlogt temp(:)];
end
logt = newlogt;

pattern = [];
for i = 1:length(blks)
    pattern = blkdiag(pattern,ones(blks(i)));
end

F_struc = [lint;logt(find(pattern),:)];
K.s(end) = [];
K.s = [K.s blks];
K.m = blks;
