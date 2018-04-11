function output = callsdpt34(interfacedata)

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

if any(K.m > 0)
    % Messy to keep track of
    options.sdpt3.smallblkdim = 0;
end

if ~isempty(interfacedata.lowrankdetails)
    options.sdpt3.smallblkdim = 1;
end

if isempty(K.schur_funs)
    % Simple...
    [blk,A,C,b,oldKs]=sedumi2sdpt3(F_struc(:,1),F_struc(:,2:end),c,K,options.sdpt3.smallblkdim);
else
    % A bit messy if we have a Schur compiler (i.e. STRUL)
    % SDPT3 reorders SDP constraints in order to put many small ones in a
    % common block. Hence, it might happen that it mixes up SDP constraints
    % with Schur compilers and those without. At the moment, we take a
    % conservative approach. If all SDP constraints have a Schur compiler,
    % we allow blocking. If not, we don't. This way our code will work    
     if ~any(cellfun(@isempty,K.schur_funs))
         % All have Schur functions
         [blk,A,C,b,oldKs]=sedumi2sdpt3(F_struc(:,1),F_struc(:,2:end),c,K,options.sdpt3.smallblkdim);        
     else
         % Messy case, don't allow blocking for now     
         options.sdpt3.smallblkdim = 1;
         [blk,A,C,b,oldKs]=sedumi2sdpt3(F_struc(:,1),F_struc(:,2:end),c,K,options.sdpt3.smallblkdim);        
     end
end
options.sdpt3.printyes=double(options.verbose);
options.sdpt3.printlevel=double(options.verbose)*3;
options.sdpt3.expon=options.sdpt3.expon(1);

% Setup the logarithmic barrier cost. We exploit the fact that we know that
% the only logaritmic cost is in the last SDP constraint
if abs(K.m) > 0
    lpLogsStart = 1;
    for i = 1:size(blk,1)
        if isequal(blk{i,1},'l')
            options.sdpt3.parbarrier{i,1} = zeros(1,blk{i,2});
        elseif isequal(blk{i,1},'u')
            lpLogsStart = 2;
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
            options.sdpt3.parbarrier{lpLogsStart,1}(end-lp_count+1) = -K.maxdetgain(i);
            lp_count = lp_count-1;
        elseif K.m(i) > 1
            % We placed it in the SDP cone
            options.sdpt3.parbarrier{end-sdp_count+1,1} = -K.maxdetgain(i);
            sdp_count = sdp_count-1;
        end
    end
    %options.saveduals = 0;
end

% Setup structures for user-defined Schur compilers
if isfield(K,'schur_funs')
    top = 1;
    if ~isempty(K.schur_funs)
        if K.f>0
            options.sdpt3.schurfun{top} = '';
            options.sdpt3.schurfun_par{top,1} = [];
            top = top+1;
        end
        if K.l > 0
            options.sdpt3.schurfun{top} = '';
            options.sdpt3.schurfun_par{top,1} = [];
            top = top+1;
        end
        if K.q > 0
            options.sdpt3.schurfun{top} = '';
            options.sdpt3.schurfun_par{top,1} = [];
            top = top+1;
        end
        if 0
            for i = 1:length(K.s)
                if ~isempty(K.schur_funs{i})
                    options.sdpt3.schurfun{top} = 'schurgateway';
                    S = createSchurFun(options,K,interfacedata,i);
                    options.sdpt3.schurfun_par{top,1} = S;
                      V = {S.extra,S.data{:}};
                      feval(S.fun,[],[],V{:});
                else
                    options.sdpt3.schurfun{top} = '';
                    options.sdpt3.schurfun_par{top,1} = [];
                end
                top = top+1;
            end
        else
            iSDP = 1;
            while top <= size(blk,1)
                for j = 1:length(blk{top,2})
                    i = oldKs(iSDP);iSDP = iSDP + 1;
                    if ~isempty(K.schur_funs{i})
                        Sname = 'schurgateway';
                     %   options.sdpt3.schurfun{top} = 'schurgateway';
                        S = createSchurFun(options,K,interfacedata,i);
                        S.j = j;
                        S.blk = blk(top,2);
                        options.sdpt3.schurfun_par{top,j} = S;
                    else
                        Sname = '';
                       % options.sdpt3.schurfun{top} = '';
                        options.sdpt3.schurfun_par{top,1} = [];
                    end
                   
                end
                options.sdpt3.schurfun{top} = Sname;                
                top = top+1;                 
            end
        end
        
    end
end

if options.savedebug
    ops = options.sdpt3;
    save sdpt3debug blk A C b ops x0 -v6
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = tic;
[obj,X,y,Z,info,runhist] =  sdpt3(blk,A,C,b,options.sdpt3,[],x0,[]);            
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

if any(K.m > 0)
   % Dual = [];
end

Primal = -y;  % Primal variable in YALMIP

% Convert error code
switch info.termcode
    case 0
        problem = 0; % No problems detected
    case {-1,-5,-9} 
        problem = 5; % Lack of progress
    case {3,-2,-3,-4,-7}
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



function S = createSchurFun(options,K,interfacedata,i);
S.extra.par = options.sdpt3;    
S.data = K.schur_data{i};
[init,loc] = ismember(K.schur_variables{i},interfacedata.used_variables);
S.index = loc;
S.fun =  K.schur_funs{i};
S.nvars = length(interfacedata.used_variables);
%options.sdpt3.schurfun_par{top,1} = S;
%  V = {S.extra,S.data{:}};
%  feval(S.fun,[],[],V{:});
