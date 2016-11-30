function output = callmosek(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
Q       = model.Q;
K       = model.K;
x0      = model.x0;
integer_variables = model.integer_variables;
binary_variables = model.binary_variables;
extended_variables = model.extended_variables;
ub      = model.ub;
lb      = model.lb;
mt      = model.monomtable;

% *********************************
% What type of variables do we have
% *********************************
model.linear_variables = find((sum(abs(mt),2)==1) & (any(mt==1,2)));
model.nonlinear_variables = setdiff((1:size(mt,1))',model.linear_variables);
model.sigmonial_variables = find(any(0>mt,2) | any(mt-fix(mt),2));

% Some meta-solver thinks we handle binaries
if ~isempty(model.binary_variables)
    integer_variables = union(model.integer_variables, model.binary_variables);
    if isempty(lb)
        lb = repmat(-inf,1,length(model.c));
    end
    lb(model.binary_variables) = max(lb(model.binary_variables),0);    
    if isempty(ub)
        ub = repmat(inf,1,length(model.c));
    end
    ub(model.binary_variables) = min(ub(model.binary_variables),1);
    model.lb = lb;
    model.ub = ub;
end

% Some meta solvers might construct model with empty cones
if any(model.K.s) && any(model.K.s == 0)
    model.K.s(model.K.s==0)=[];
end
if any(model.K.q) && any(model.K.q == 0)
    model.K.q(model.K.q==0)=[];
end

if ~isempty(model.sigmonial_variables) | isequal(model.solver.version,'GEOMETRIC')
    [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_geometric(model);          
else
    [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqpsocpsdp(model);        
end

infostr = yalmiperror(problem,'MOSEK');	

% Save all data sent to solver?
if options.savesolverinput
    solverinput.prob = prob;
    solverinput.param = options.mosek;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.r = r;
    solveroutput.res = res;
    
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);


function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqpsocpsdp(model);

if nnz(model.Q)==0 & isempty(model.integer_variables) && isempty(model.x0)
    % Standard cone problem
    [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqpsocpsdpdual(model);
    return
elseif model.K.s(1)>0
    % Semidefinite program with integer variables or quadratic objective
    % or specified initial point
    [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_sdp(model);
else
    % Integer conic program + possibly quadratic objective
    [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqp(model);        
end


function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_geometric(model);
   
x = [];
D_struc = [];
r = [];
res = [];
solvertime = 0;
[prob,problem] = yalmip2geometric(model.options,model.F_struc,model.c,model.Q,model.K,model.ub,model.lb,model.monomtable,model.linear_variables,model.extended_variables);
if problem == 0

    % Mosek does not support equalities
    if ~isempty(prob.G)
        prob.A = [prob.A;prob.G;-prob.G];
        prob.b = [prob.b;prob.h;1./prob.h];
        prob.map = [prob.map;max(prob.map) + (1:2*length(prob.h))'];
    end
    if model.options.savedebug
        save mosekdebug prob
    end

    param = model.options.mosek;

    % Call MOSEK   
    showprogress('Calling MOSEK',model.options.showprogress);    
    if model.options.verbose == 0  
        solvertime = tic;
        res = mskgpopt(prob.b,prob.A,prob.map,param,'minimize echo(0)');
        solvertime = toc(solvertime);
    else
        solvertime = tic;
    	res = mskgpopt(prob.b,prob.A,prob.map,param,'minimize');     
        solvertime = toc(solvertime);
    end
    sol = res.sol;
    
    x = zeros(length(model.c),1);
    x(model.linear_variables)=exp(res.sol.itr.xx);
    D_struc = [];

    % Check, currently not exhaustive...
    switch res.sol.itr.prosta
        case 'PRIMAL_AND_DUAL_FEASIBLE'
            problem = 0;
        case 'PRIMAL_INFEASIBLE'
            problem = 1;
        case 'DUAL_INFEASIBLE'
            problem = 2;
        case 'UNKNOWN'
            problem = 9;
        otherwise
            problem = -1;
    end
end

function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqpsocpsdpdual(model)

% Convert if the caller is bnb or bmibnb which might have appended bounds
% Sure, we could setup model with bounds, but... 
[model.F_struc,model.K] = addStructureBounds(model.F_struc,model.K,model.ub,model.lb);

param = model.options.mosek;

prob.c = model.F_struc(1:model.K.f+model.K.l+sum(model.K.q),1);
prob.a = -model.F_struc(1:model.K.f+model.K.l+sum(model.K.q),2:end)';
prob.blc = -model.c;
prob.buc = -model.c;
prob.blx = -inf(size(prob.a,2),1);
prob.bux = inf(size(prob.a,2),1);
top = model.K.f+model.K.l;
prob.blx(1+model.K.f:model.K.f+model.K.l) = 0;

if model.K.q(1)>0
    for i = 1:length(model.K.q)
        prob.cones{i}.type = 'MSK_CT_QUAD';
        prob.cones{i}.sub  = top+1:top+model.K.q(i);
        top = top + model.K.q(i);
    end
end

if model.K.s(1)>0
   prob = appendSDPdata(model.F_struc,model.K,prob);
end

if model.options.savedebug
    ops = model.options;
    save mosekdebug prob param ops
end

[r,res,solvertime] = doCall(prob,param,model.options);

try
    x = res.sol.itr.y;
catch
    x = nan(length(model.c),1);    
end

if model.options.saveduals & ~isempty(x)
    try       
        D_struc_SDP = zeros(sum(model.K.s.^2),1);
        top = 1;
        dtop = 1;
        for i = 1:length(model.K.s)          
            n = model.K.s(i);
            I = find(tril(ones(n)));
            v = res.sol.itr.barx(top:((top+n*(n+1)/2)-1));
            D_struc_SDP(dtop + I - 1) = v;
            in = ceil(I/n);
            jn = mod(I-1,n)+1;
            D_struc_SDP(dtop + (jn-1)*n+in - 1) = v;
            top = top + n*(n+1)/2;
            dtop = dtop  + n^2;
        end
        D_struc = [res.sol.itr.xx;D_struc_SDP];
    catch
        D_struc = [];
    end
else
    D_struc = [];
end

problem = MosekYALMIPError(res);

function [res,sol,solvertime] = doCall(prob,param,options)

showprogress('Calling Mosek',options.showprogress);
if options.verbose == 0
    solvertime = tic;
    [res,sol] = mosekopt('minimize echo(0)',prob,param);    
    solvertime = toc(solvertime);
else
    solvertime = tic;
    [res,sol] = mosekopt('minimize info',prob,param);
    solvertime = toc(solvertime);
end

function problem = MosekYALMIPError(res)

if res.rcode == 2001
    problem = 1;
    return
elseif res.rcode == 10007
    problem = 16;
    return
end

switch res.sol.itr.prosta
    case 'PRIMAL_AND_DUAL_FEASIBLE'        
        problem = 0;
    case 'DUAL_INFEASIBLE'
        problem = 1;
    case 'PRIMAL_INFEASIBLE'
        problem = 2;
    case 'MSK_RES_TRM_USER_CALLBACK'
        problem = 16;
    case 'MSK_RES_TRM_STALL'
        problem = 4;
    case 'UNKNOWN'
        try
            if isequal(res.rcodestr,'MSK_RES_TRM_STALL')
                problem = 4;
            else
                problem = 9;
            end
        catch
            problem = 9;
        end
    otherwise
        problem = -1;
end

function model = appendBounds(model);

FstrucNew = [];
if ~isempty(model.lb)
    ii = find(~isinf(model.lb));
    FstrucNew = [-model.lb(ii) sparse(1:length(ii),ii,size(model.F_struc,2)-1,length(ii))];
end
if ~isempty(model.ub)
    ii = find(~isinf(model.ub));
    FstrucNew = [model.lb(ii) -sparse(1:length(ii),ii,size(model.F_struc,2)-1,length(ii))];
end

function prob = appendSDPdata(F_struc,K,prob)

prob.bardim = K.s;
prob.barc.subj = [];
prob.barc.subk = [];
prob.barc.subl = [];
prob.barc.val = [];
prob.bara.subi = [];
prob.bara.subj = [];
prob.bara.subk = [];
prob.bara.subl = [];
prob.bara.val = [];

C = F_struc(:,1);
A = -F_struc(:,2:end);

% -- Faster fix by Shahar
tops = [1 cumsum(K.s.^2)+1];
top = 1+K.f+K.l+sum(K.q);
[ii,jj,kk] = find(A(top:top + sum(K.s.^2)-1,:));
allcon = floor(interp1(tops,1:length(tops),ii,'linear'));
all_iilocal = ii-tops(allcon)'+1;
a = all_iilocal;
b = K.s(allcon);
allcol = ceil(a(:)./b(:))';
allrow = a(:)' - (allcol-1).*b(:)';
allvar = jj;
allval = kk;
% sort (for backward compatibility?)
[~,ind_sort] = sort(allcon);
allcon = allcon(ind_sort)';
allcol = allcol(ind_sort)';
allrow = allrow(ind_sort)';
allvar = allvar(ind_sort)';
allval = allval(ind_sort)';
% --

keep = find(allrow >= allcol);
allcol = allcol(keep);
allrow = allrow(keep);
allcon = allcon(keep);
allvar = allvar(keep);
allval = allval(keep);
prob.bara.subi = [prob.bara.subi allvar];
prob.bara.subj = [prob.bara.subj allcon];
prob.bara.subk = [prob.bara.subk allrow];
prob.bara.subl = [prob.bara.subl allcol];
prob.bara.val = [prob.bara.val allval];

for j = 1:length(K.s)
    n = K.s(j);
    Ci = C(top:top+n^2-1);
    Ci = tril(reshape(Ci,n,n));
    [k,l,val] = find(Ci);
    prob.barc.subj = [prob.barc.subj j*ones(1,length(k))];
    prob.barc.subk = [prob.barc.subk k(:)'];
    prob.barc.subl = [prob.barc.subl l(:)'];
    prob.barc.val = [prob.barc.val val(:)'];     
    top = top + n^2;
end

function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_sdp(model)

param = model.options.mosek;

% Convert
[model.F_struc,model.K] = addStructureBounds(model.F_struc,model.K,model.ub,model.lb);
prob = yalmip2SDPmosek(model);

% Debug?
if model.options.savedebug
    save mosekdebug prob param
end

if model.options.verbose == 0
    solvertime = tic;
    [r,res] = mosekopt('minimize echo(0)',prob,param);    
    solvertime = toc(solvertime);
else
    solvertime = tic;
    [r,res] = mosekopt('minimize info',prob,param);
    solvertime = toc(solvertime);
end

if res.rcode == 2001
    res.sol.itr.prosta = 'DUAL_INFEASIBLE';
end

try
    x = res.sol.itr.y;
catch
    x = nan(length(model.c),1);    
end

if model.options.saveduals & ~isempty(x)
    D_struc = [res.sol.itr.xx];   
    D_struc_SDP = sum(model.K.s.^2);
    top = 1;
    dtop = 1;
    for i = 1:length(model.K.s)
        X = zeros(model.K.s(i));
        n = model.K.s(i);
        I = find(tril(ones(n)));
        D_struc_SDP(dtop + I - 1) = res.sol.itr.barx(top:((top+n*(n+1)/2)-1));
        X(I) = res.sol.itr.barx(top:((top+n*(n+1)/2)-1));
        X = X + tril(X,-1)';
        D_struc = [D_struc;X(:)];
        top = top + n*(n+1)/2;
        dtop = dtop  + n^2;
    end
else
    D_struc = [];
end

switch res.sol.itr.prosta
    case 'PRIMAL_AND_DUAL_FEASIBLE'
        problem = 0;
    case 'DUAL_INFEASIBLE'
        problem = 1;
    case 'PRIMAL_INFEASIBLE'
        problem = 2;
    case 'UNKNOWN'
        try
            if isequal(res.rcodestr,'MSK_RES_TRM_STALL')
                problem = 4;
            else
                problem = 9;
            end
        catch
            problem = 9;
        end
    otherwise
        problem = -1;
end

function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqp(model);

prob.c = model.c;
if ~isempty(model.F_struc)
    prob.a = -model.F_struc(:,2:end);
    prob.buc = full(model.F_struc(:,1));
    prob.blc = repmat(-inf,length(prob.buc),1);
else
    prob.a = sparse(ones(1,length(model.c))); % Dummy constraint
    prob.buc = inf;
    prob.blc = -inf;
end
if isempty(model.lb)
    prob.blx = repmat(-inf,1,length(model.c));
else
    prob.blx = model.lb;
end
if isempty(model.ub)
    prob.bux = repmat(inf,1,length(model.c));
else
    prob.bux = model.ub;
end

if model.K.f>0
    prob.blc(1:model.K.f) = prob.buc(1:model.K.f);
end

[prob.qosubi,prob.qosubj,prob.qoval] = find(tril(sparse(2*model.Q)));    

if model.K.q(1)>0
    nof_new = sum(model.K.q);
    prob.a = [prob.a [spalloc(model.K.f,nof_new,0);spalloc(model.K.l,nof_new,0);speye(nof_new);spalloc(sum(model.K.s.^2),nof_new,0)]];
    prob.blc(1+model.K.f+model.K.l:model.K.f+model.K.l+sum(model.K.q)) = prob.buc(1+model.K.f+model.K.l:model.K.f+model.K.l+sum(model.K.q)); % Note, fixed the SDP cones too
      
    prob.c = [prob.c;zeros(nof_new,1)];
    top = size(model.F_struc,2)-1;
    for i = 1:length(model.K.q)
        prob.cones{i}.type = 'MSK_CT_QUAD';
        prob.cones{i}.sub  = top+1:top+model.K.q(i);
        prob.blx(top+1:top+model.K.q(i)) = -inf;
        prob.bux(top+1:top+model.K.q(i)) = inf;        
        top = top + model.K.q(i);
    end    
end

if ~isempty(model.integer_variables)
    prob.ints.sub = model.integer_variables;
end

param = model.options.mosek;

if ~isempty(model.x0)
    if model.options.usex0
        prob.sol.int.xx = zeros(max([length(model.Q) size(prob.a,2)]),1);
        prob.sol.int.xx(model.integer_variables) = model.x0(model.integer_variables);
        evalc('[r,res] = mosekopt (''symbcon'')');
        sc = res.symbcon ;
        param.MSK_IPAR_MIO_CONSTRUCT_SOL = sc.MSK_ON;
    end
end

% Debug?
if model.options.savedebug
    save mosekdebug prob param
end

% Call MOSEK
showprogress('Calling MOSEK',model.options.showprogress);
if model.options.verbose == 0
    solvertime = tic;
    [r,res] = mosekopt('minimize echo(0)',prob,param); 
    solvertime = toc(solvertime);
else
    solvertime = tic;
    [r,res] = mosekopt('minimize',prob,param);     
    solvertime = toc(solvertime);
end

if (r == 1010) || (r == 1011) | (r==1001)
    problem = -5;
    x = [];
    D_struc = [];
elseif r == 1200
    problem = 7;
    x = [];
    D_struc = [];
else
    % Recover solutions
    sol = res.sol;
    if isempty(model.integer_variables)
        x = sol.itr.xx(1:length(model.c)); % Might have added new ones                
        D_struc = (sol.itr.suc-sol.itr.slc);        
        error_message = sol.itr.prosta;
    else        
        try
        x = sol.int.xx(1:length(model.c)); % Might have added new ones
        D_struc = [];
        error_message = sol.int.prosta;
        catch
            x = [];
            error_message = 'crash';
            D_struc = [];
        end
      
    end
    
    switch error_message
    case {'PRIMAL_AND_DUAL_FEASIBLE','PRIMAL_FEASIBLE'}
        problem = 0;
    case 'PRIMAL_INFEASIBLE'
        problem = 1;
    case 'DUAL_INFEASIBLE'        
        problem = 2;
    case 'PRIMAL_INFEASIBLE_OR_UNBOUNDED'
        problem = 12;
    otherwise
        problem = -1;
    end
end
