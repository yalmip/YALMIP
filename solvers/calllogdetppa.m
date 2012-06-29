function output = calllogdetppa(interfacedata)

% Author Johan Löfberg
% $Id: calllogdetppa.m,v 1.21 2010-01-13 13:49:21 joloef Exp $ 

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

% SDPNAL does not support multiple blocks
options.sdpt3.smallblkdim = 0;

% Convert from internal (sedumi-like) format
[blk,A,C,b,oldKs]=sedumi2sdpt3(F_struc(:,1),F_struc(:,2:end),c,K,options.sdpt3.smallblkdim);

% Setup the logarithmic barrier cost. We exploit the fact that we know that
% the only logaritmic cost is in the last SDP constraint
if abs(K.m) > 0
    for i = 1:size(blk,1)
        if isequal(blk{i,1},'l')
            options.sdpt3.parbarrier{i,1} = zeros(1,blk{i,2});
        else
            options.sdpt3.parbarrier{i,1} = 0*blk{i,2};
        end
    end
    n_sdp_logs = nnz(K.m > 1);
    n_lp_logs  = nnz(K.m == 1);
    if n_lp_logs>0
        lp_count = n_lp_logs;
    end
    if n_sdp_logs>0
        sdp_count = n_sdp_logs;
    end  
    for i = 1:length(K.m)
        if K.m(i) == 1
            % We placed it in the linear cone
            options.sdpt3.parbarrier{1,1}(end-lp_count+1) = -K.maxdetgain(i);
            lp_count = lp_count-1;
        elseif K.m(i) > 1
            % We placed it in the SDP cone
            options.sdpt3.parbarrier{end-sdp_count+1,1} = -K.maxdetgain(i);
            sdp_count = sdp_count-1;
        end
    end
    %options.saveduals = 0;
   mu0 = [options.sdpt3.parbarrier{:}];
else
   mu0 = zeros(K.l+length(find(K.s)),1);
end

%options.logdetppa.mu0 = [options.sdpt3.parbarrier{:}];
if options.savedebug
    ops = options.logdetppa;
    save logdetppadebug blk A C b ops x0 -v6
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = clock;
if options.verbose==0 % SDPT3 does not run silent despite printyes=0!
   evalc('[obj,X,y,Z,info,runhist] =   logdetPPA(blk,A,C,b,mu0,options.logdetppa);');
else    
    [obj,X,y,Z,info,runhist] =  logdetPPA(blk,A,C,b,mu0,options.logdetppa);            
end

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

if any(K.m > 0)
   % Dual = [];
end

% if options.removethem
% Dual = [];
% end

solvertime = etime(clock,solvertime);
Primal = -y;  % Primal variable in YALMIP
problem = -7;
% 
% % Convert error code
% switch info.termcode
%     case -1
%         problem = 4;
%     case 0
%         problem = 0;
%     case 1
%         problem = 5;
%     case 2
%         problem = 3;   
%     otherwise
%         problem = 11;
% end

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
output.Primal      = Primal;
output.Dual        = Dual;
output.Slack       = Slack;
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;

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
