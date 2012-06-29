function output = callpenbmi(interfacedata);

% Author Johan Löfberg
% $Id: callpenbmim.m,v 1.24 2007-07-28 14:09:01 joloef Exp $

if any(interfacedata.variabletype > 2)
    % Polynomial problem, YALMIP has to bilienarize
    interfacedata.high_monom_model=[];
    output = callpenbmi_with_bilinearization(interfacedata);
else
    % Old standard code
    output = callpenbmi_without_bilinearization(interfacedata);
end

function output = callpenbmi_without_bilinearization(interfacedata);

% Retrieve needed data
clear penbmim

options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
monomtable = interfacedata.monomtable;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

% Linear before bilinearize
temp = sum(monomtable,2)>1;
tempnonlinearindicies = find(temp);
templinearindicies = find(~temp);

% Any stupid constant>0 constraints
% FIX : Recover duals afterwards
% Better fix : Do this outside
zrow = [];
if K.l >0
    zrow = find(any(F_struc(1:K.l+K.f,:),2)==0);
    if ~isempty(zrow)
        K.l = K.l - nnz(zrow>K.f);
        K.f = K.f - nnz(zrow<=K.f);       
        F_struc(zrow,:) = [];
    end
end

% Bounded variables converted to constraints
if ~isempty(ub)
    [F_struc,K] = addbounds(F_struc,K,ub,lb);
end

% This one only occurs if called from bmibnb
if K.f>0
    F_struc = [-F_struc(1:K.f,:);F_struc];
    F_struc(1:K.f,1) = F_struc(1:K.f,1)+sqrt(eps);
    K.l = K.l + 2*K.f;
    K.f = 0;
end

if isempty(monomtable)
    monomtable = eye(length(c));
end

temp = sum(monomtable,2)>1;
nonlinearindicies = find(temp);
linearindicies = find(~temp);
c0 = c;
c  = c(linearindicies);
Q = Q(linearindicies,linearindicies);
nonlinear_scalars = [];

% Any non-linear scalar inequalities?
% Move these to the BMI part
if K.l>0
    nonlinear_scalars = find(any(full(F_struc(1:K.l,[nonlinearindicies(:)'+1])),2));
    if ~isempty(nonlinear_scalars)
        Kold = K;
        % SETDIFF DONE FASTER
         aa = 1:K.l;
         bb = nonlinear_scalars;
         tf = ~(ismembc(aa,bb)); 
         cc = aa(tf); 
         cc = unique(cc); 
         linear_scalars = cc;
%        linear_scalars = setdiff1(1:K.l,nonlinear_scalars);
        F_struc = [F_struc(linear_scalars,:);F_struc(nonlinear_scalars,:);F_struc(K.l+1:end,:)];     
        K.l = K.l-length(nonlinear_scalars);
        if (length(K.s)==1) & (K.s==0)
            K.s    = [repmat(1,1,length(nonlinear_scalars))];
            K.rank = repmat(1,1,length(nonlinear_scalars)); 
        else
            K.s    = [repmat(1,1,length(nonlinear_scalars)) K.s];
            K.rank = [repmat(1,1,length(nonlinear_scalars)) K.rank];
        end
    end
end

if ~isempty(F_struc)
    pen = sedumi2pen(F_struc(:,[1 linearindicies(:)'+1]),K,c,x0);
else
    pen = sedumi2pen([],K,c,x0);
end

if ~isempty(nonlinearindicies)
    bmi = sedumi2pen(F_struc(:,[nonlinearindicies(:)'+1]),K,[],[]);
    pen.ki_dim = bmi.ai_dim;
    % Nonlinear index
    pen.ki_dim = bmi.ai_dim;
    pen.ki_row = bmi.ai_row;
    pen.ki_col = bmi.ai_col;
    pen.ki_nzs = bmi.ai_nzs;
    pen.ki_val = bmi.ai_val;
    if 0
        for i = 1:length(bmi.ai_idx)
            nl = nonlinearindicies(1+bmi.ai_idx(i));
            v = find(monomtable(nl,:));
            if length(v)==1
                v(2)=v(1);
            end
            pen.ki_idx(i)=find(linearindicies == v(1));
            pen.kj_idx(i)=find(linearindicies == v(2));
        end
    else
        top = 1;
        [ii,jj,kk] = find(monomtable(nonlinearindicies(1+bmi.ai_idx),:)');
        pen.ki_idx = zeros(1,length(bmi.ai_idx));
        pen.kj_idx = zeros(1,length(bmi.ai_idx));
        for i = 1:length(bmi.ai_idx)
            if kk(top)==2
                v(1) = ii(top);
                v(2) = ii(top);
                top = top + 1;
            else
                v(1) = ii(top);
                v(2) = ii(top+1);
                top = top + 2;
            end
            % FIX : precompute this map
            pen.ki_idx(i)=find(linearindicies == v(1));
            pen.kj_idx(i)=find(linearindicies == v(2));
        end
    end

else
    pen.ki_dim = 0*pen.ai_dim;
    pen.ki_row = 0;
    pen.ki_col = 0;
    pen.ki_nzs = 0;
    pen.ki_idx = 0;
    pen.kj_idx = 0;
    pen.kj_val = 0;    
end

if nnz(Q)>0
    [row,col,vals] = find(triu(Q));
    pen.q_nzs = length(row);
    pen.q_val = vals';
    pen.q_col = col'-1;
    pen.q_row = row'-1;
else
    pen.q_nzs = 0;
    pen.q_val = 0;
    pen.q_col = 0;
    pen.q_row = 0;    
end

ops = struct2cell(options.penbmi);ops = [ops{1:end}];
pen.ioptions = ops(1:12);
pen.foptions = ops(13:end);
pen.ioptions(4) = max(0,min(3,options.verbose+1));
if pen.ioptions(4)==1
    pen.ioptions(4)=0;
end

% ****************************************
% UNCOMMENT THIS FOR PENBMI version 1
% ****************************************
%pen.ioptions = pen.ioptions(1:8);
%pen.foptions = pen.foptions(1:8);

if ~isempty(x0)    
    pen.x0(isnan(pen.x0)) = 0;
    pen.x0 = x0(linearindicies);    
    pen.x0 = pen.x0(:)';
end

if options.savedebug
    save penbmimdebug pen
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end

solvertime = clock; 
[f,xout,u,iflag,niter,feas] = penbmim(pen);
solvertime = etime(clock,solvertime);

if options.saveduals & isempty(zrow)
    
    % Get dual variable
    % First, get the nonlinear scalars treated as BMIs
    if ~isempty(nonlinear_scalars)
        n_orig_scalars = length(nonlinear_scalars)+K.l;
        linear_scalars = setdiff(1:n_orig_scalars,nonlinear_scalars);
        u_nonlinear=u(K.l+1:K.l+length(nonlinear_scalars));
        u(K.l+1:K.l+length(nonlinear_scalars))=[];
        u_linear = u(1:K.l);
        u_scalar = zeros(1,n_orig_scalars);
        u_scalar(linear_scalars)=u_linear;
        u_scalar(nonlinear_scalars)=u_nonlinear;
        u = [u_scalar u(1+K.l:end)];
        K = Kold;
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
else
    D_struc = [];
end

%Recover solution
if isempty(nonlinearindicies)
    x = xout(:);
else
    x = zeros(length(interfacedata.c),1);    
    for i = 1:length(templinearindicies)
        x(templinearindicies(i)) = xout(i);
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
infostr = yalmiperror(problem,interfacedata.solver.tag);	

if options.savesolveroutput   
    solveroutput.f = f;
    solveroutput.xout = xout;
    solveroutput.u = u;
    solveroutput.iflag = iflag;
    solveroutput.niter = niter;
    solveroutput.feas = feas;
else
    solveroutput = [];
end

if options.savesolverinput
    solverinput.pen = pen;
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

function output = callpenbmi_with_bilinearization(interfacedata);

% Bilinearize
[p,changed] = convert_polynomial_to_quadratic(interfacedata);

% Convert bilinearizing equalities to inequalities
if p.K.f>0
    p.F_struc = [-p.F_struc(1:p.K.f,:);p.F_struc];
    p.F_struc(1:p.K.f,1) = p.F_struc(1:p.K.f,1)+sqrt(eps);
    p.K.l = p.K.l + 2*p.K.f;
    p.K.f = 0;
end

% Solve bilinearized problem
 output = callpenbmi_without_bilinearization(p);
 
 % Get our original variables & duals
 output.Primal = output.Primal(1:length(interfacedata.c));
 if ~isempty(output.Dual)
     n_equ = p.K.f - interfacedata.K.f;
     % First 2*n_eq are the duals for the new inequalities
     output.Dual = output.Dual(1+2*n_equ:end);
 end
 
 
 


