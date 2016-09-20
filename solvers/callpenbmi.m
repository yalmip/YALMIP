function output = callpenbmi(interfacedata);

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
monomtable = interfacedata.monomtable;
ub      = interfacedata.ub;
lb      = interfacedata.lb;


% Bounded variables converted to constraints
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

if K.f>0
    F_struc = [-F_struc(1:K.f,:);F_struc];
    F_struc(1:K.f,1) = F_struc(1:K.f,1)+sqrt(eps);
    K.l = K.l + 2*K.f;
    K.f = 0;
end

nonlinearindicies = find(sum(monomtable,2)>1);
linearindicies = setdiff(1:length(c),nonlinearindicies);
c0 = c;
c  = c(linearindicies);

% Any non-linear scalar inequalities?
% Move these to the BMI part
if K.l>0
    nonlinear_scalars = find(any(full(F_struc(1:K.l,[nonlinearindicies(:)'+1])),2));
    if ~isempty(nonlinear_scalars)
        Kold = K;
        linear_scalars = setdiff(1:K.l,nonlinear_scalars);
        F_struc = [F_struc(linear_scalars,:);F_struc(nonlinear_scalars,:);F_struc(K.l+1:end,:)];
        K.l = K.l-length(nonlinear_scalars);
        if (length(K.s)==1) & (K.s==0)
            K.s = [repmat(1,1,length(nonlinear_scalars))];
        else
            K.s = [repmat(1,1,length(nonlinear_scalars)) K.s];
        end
    end
end

if ~isempty(F_struc)
    penstruct = sedumi2pen(F_struc(:,[1 linearindicies(:)'+1]),K,c,x0);
else
    penstruct = sedumi2pen([],K,c,x0);
end

if ~isempty(nonlinearindicies)
    bmi = sedumi2pen(F_struc(:,[nonlinearindicies(:)'+1]),K,[],[]);
    penstruct.ki_dim = bmi.ai_dim;
    % Nonlinear index
    penstruct.ki_dim = bmi.ai_dim;
    penstruct.ki_row = bmi.ai_row;
    penstruct.ki_col = bmi.ai_col;
    penstruct.ki_nzs = bmi.ai_nzs;
    penstruct.ki_val = bmi.ai_val;
    for i = 1:length(bmi.ai_idx)
        nl = nonlinearindicies(1+bmi.ai_idx(i));
        v = find(monomtable(nl,:));
        if length(v)==1
            v(2)=v(1);
        end

        penstruct.ki_idx(i)=v(1);
        penstruct.kj_idx(i)=v(2);

    end
else
    penstruct.ki_dim = 0*penstruct.ai_dim;
    penstruct.ki_row = [];
    penstruct.ki_col = [];
    penstruct.ki_nzs = [];
    penstruct.ki_val = [];
    penstruct.ki_idx = [];
    penstruct.kj_idx = [];
    penstruct.kj_val = [];
end

if nnz(Q)>0
    [row,col,vals] = find(triu(Q));
    penstruct.q_nzs = length(row);
    penstruct.q_val = vals;
    penstruct.q_col = col-1;
    penstruct.q_row = row-1;
else
    penstruct.q_nzs = 0;
    penstruct.q_val = 0;
    penstruct.q_col = 0;
    penstruct.q_row = 0;
end

ops = struct2cell(options.penbmi);ops = [ops{1:end}];
penstruct.ioptions = ops(1:12);
penstruct.foptions = ops(13:end);
penstruct.ioptions(4) = max(0,min(3,options.verbose+1));
if penstruct.ioptions(4)==1
    penstruct.ioptions(4)=0;
end

% ****************************************
% UNCOMMENT THIS IF USING PENBMI version 1
% ****************************************
% penstruct.ioptions = penstruct.ioptions(1:8);
% penstruct.foptions = penstruct.foptions(1:8);

if ~isempty(x0)
    penstruct.x0 = x0(linearindicies);
    penstruct.x0 = penstruct.x0(:)';
end

% FIX
if penstruct.mconstr == 0
    penstruct.msizes = [];
end

if options.savedebug
    save penbmidebug penstruct
end

showprogress('Calling PENBMI',options.showprogress);
solvertime = tic;
try    
    if all(c==0)
        [xout, fx, u, iresults, fresults, iflag] = pen(penstruct,1);
    else
        [xout, fx, u, iresults, fresults, iflag] = pen(penstruct,0);
    end
catch
    % Fix for bug i tomlab
    if all(c==0)
        [xout, fx, u, iresults, fresults, iflag] = pen(penstruct);
    else
        [xout, fx, u, iresults, fresults, iflag] = pen(penstruct);
    end    
end
solvertime = toc(solvertime);

% Get dual variable
% First, get the nonlinear scalars treated as BMIs
if exist('nonlinear_scalars')
    if ~isempty(nonlinear_scalars)
        u = u(:);
        n_orig_scalars = length(nonlinear_scalars)+K.l;
        linear_scalars = setdiff(1:n_orig_scalars,nonlinear_scalars);
        u_nonlinear=u(K.l+1:K.l+length(nonlinear_scalars));
        u(K.l+1:K.l+length(nonlinear_scalars))=[];
        u_linear = u(1:K.l);
        u_scalar = zeros(1,n_orig_scalars);
        u_scalar(linear_scalars)=u_linear;
        u_scalar(nonlinear_scalars)=u_nonlinear;
        u = [u_scalar(:);u(1+K.l:end)];
        K = Kold;
    end
end

u = u(:);
D_struc = u(1:1:K.l);
if length(K.s)>0
    if K.s(1)>0
        pos = K.l+1;
        for i = 1:length(K.s)
            temp = zeros(K.s(i),K.s(i));
            vecZ = u(pos:pos+0.5*K.s(i)*(K.s(i)+1)-1);
            top = 1;
            for j = 1:K.s(i)
                len = K.s(i)-j+1;
                temp(j:end,j)=vecZ(top:top+len-1);top=top+len;
            end
            temp = (temp+temp');j = find(speye(K.s(i)));temp(j)=temp(j)/2;
            D_struc = [D_struc;temp(:)];
            pos = pos + (K.s(i)+1)*K.s(i)/2;
        end
    end
end


%Recover solution
if isempty(nonlinearindicies)
    x = xout(:);
else
    x = zeros(length(c0),1);
    for i = 1:length(linearindicies)
        x(linearindicies(i)) = xout(i);
    end
end

problem = 0;
switch iflag
    case 0
        problem = 0; % OK
    case {1,3}
        problem = 4;
    case 2
        problem = 1; % INFEASIBLE
    case 4
        problem = 3; % Numerics
    case 5
        problem = 7;
    case {6,7}
        problem = 11;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,'PENBMI/TOMLAB');

if options.savesolveroutput
    solveroutput.xout = xout;
    solveroutput.fx = fx;
    solveroutput.u = u;
    solveroutput.iresults = iresults;
    solveroutput.fresults = fresults;
    solveroutput.iflag = iflag;
else
    solveroutput = [];
end

if options.savesolverinput
    solverinput.penstruct = penstruct;
else
    solverinput = [];
end

% Standard interface 
output = createOutputStructure(x(:),D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);