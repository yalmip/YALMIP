function output = callmosek(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
extended_variables = interfacedata.extended_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
mt      = interfacedata.monomtable;

% *********************************
% What type of variables do we have
% *********************************
linear_variables = find((sum(abs(mt),2)==1) & (any(mt==1,2)));
nonlinear_variables = setdiff((1:size(mt,1))',linear_variables);
sigmonial_variables = find(any(0>mt,2) | any(mt-fix(mt),2));

% Some meta-solver thinks we handle binaries
if ~isempty(interfacedata.binary_variables)
    integer_variables = union( interfacedata.integer_variables, interfacedata.binary_variables);
    if isempty(lb)
        lb = repmat(-inf,1,length(interfacedata.c));
    end
    lb(interfacedata.binary_variables) = max(lb(interfacedata.binary_variables),0);    
    if isempty(ub)
        ub = repmat(inf,1,length(interfacedata.c));
    end
    ub(interfacedata.binary_variables) = min(ub(interfacedata.binary_variables),1);
end

if ~isempty(sigmonial_variables) | isequal(interfacedata.solver.version,'GEOMETRIC')
    [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_geometric(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,extended_variables);          
else
    [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqpsocpsdp(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,integer_variables,x0);        
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
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;

function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqpsocpsdp(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,integer_variables,x0);

if nnz(Q)==0 & isempty(integer_variables) && isempty(x0)
    [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqpsocpsdpdual(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,integer_variables,x0);
    return
end

% SETUP a naive model interpreting directly in primal space. Will not be
% efficient for large-scale problems, but integrality can only be appended
% in this space. quadratic terms are easier here as well

prob.c = c;
if ~isempty(F_struc)
    prob.a = -F_struc(:,2:end);
    prob.buc = full(F_struc(:,1));
    prob.blc = repmat(-inf,length(prob.buc),1);
else
    prob.a = sparse(ones(1,length(c))); % Dummy constraint
    prob.buc = inf;
    prob.blc = -inf;
end
if isempty(lb)
    prob.blx = repmat(-inf,1,length(c));
else
    prob.blx = lb;
end
if isempty(ub)
    prob.bux = repmat(inf,1,length(c));
else
    prob.bux = ub;
end

if K.f>0
    prob.blc(1:K.f) = prob.buc(1:K.f);
end

[prob.qosubi,prob.qosubj,prob.qoval] = find(tril(sparse(2*Q)));    

if K.q(1)>0
    nof_new = sum(K.q);
    prob.a = [prob.a [spalloc(K.f,nof_new,0);spalloc(K.l,nof_new,0);speye(nof_new);spalloc(sum(K.s.^2),nof_new,0)]];
    prob.blc(1+K.f+K.l:K.f+K.l+sum(K.q)) = prob.buc(1+K.f+K.l:K.f+K.l+sum(K.q));    
    prob.c = [prob.c;zeros(nof_new,1)];
    
    top = size(F_struc,2)-1;
    for i = 1:length(K.q)
        prob.cones{i}.type = 'MSK_CT_QUAD';
        prob.cones{i}.sub  = top+1:top+K.q(i);
        prob.blx(top+1:top+K.q(i)) = -inf;
        prob.bux(top+1:top+K.q(i)) = inf;        
        top = top + K.q(i);
    end    
end

if ~isempty(integer_variables)
    prob.ints.sub = integer_variables;
end

param = options.mosek;

if ~isempty(x0)
    if options.usex0
        prob.sol.int.xx = zeros(max([length(Q) size(prob.a,2)]),1);
        prob.sol.int.xx(integer_variables) = x0(integer_variables);
        evalc('[res,sol] = mosekopt (''symbcon'')');
        sc = sol.symbcon ;
        param.MSK_IPAR_MIO_CONSTRUCT_SOL = sc.MSK_ON;
    end
end

% Debug?
if options.savedebug
    save mosekdebug prob param
end

% Call MOSEK
[r,res,solvertime] = doCall(prob,param,options);

x = [];
D_struc = [];
if (r == 1010) || (r == 1011) | (r==1001)
    problem = -5;
elseif r == 1200
    problem = 7;
elseif r == 10007
    problem = 16;    
else
    % Recover solutions
    sol = res.sol;
    if isempty(integer_variables)
        x = sol.itr.xx(1:length(c)); % Might have added new ones                
        D_struc = (sol.itr.suc-sol.itr.slc);        
        error_message = sol.itr.prosta;
    else        
        try
        x = sol.int.xx(1:length(c)); % Might have added new ones
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

function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_geometric(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,extended_variables);
   
x = [];
D_struc = [];
r = [];
res = [];
solvertime = 0;
[prob,problem] = yalmip2geometric(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,extended_variables);
if problem == 0

    % Mosek does not support equalities
    if ~isempty(prob.G)
        prob.A = [prob.A;prob.G;-prob.G];
        prob.b = [prob.b;prob.h;1./prob.h];
        prob.map = [prob.map;max(prob.map) + (1:2*length(prob.h))'];
    end
    if options.savedebug
        save mosekdebug prob
    end

    param = options.mosek;

    % Call MOSEK
    solvertime = clock;
    showprogress('Calling MOSEK',options.showprogress);
    solvertime = clock;
    if options.verbose == 0  
        res = mskgpopt(prob.b,prob.A,prob.map,param,'minimize echo(0)');
    else
    	res = mskgpopt(prob.b,prob.A,prob.map,param,'minimize');     
    end
    sol = res.sol;
    solvertime = etime(clock,solvertime);

    x = zeros(length(c),1);
    x(linear_variables)=exp(res.sol.itr.xx);
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

function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqpsocpsdpdual(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,integer_variables,x0)

% Convert if the caller is bnb or bmibnb which might have appended bounds
% Sure, we could setup model with bounds, but... 
[F_struc,K] = addbounds(F_struc,K,ub,lb);

param = options.mosek;

prob.c = F_struc(1:K.f+K.l+sum(K.q),1);
prob.a = -F_struc(1:K.f+K.l+sum(K.q),2:end)';
prob.blc = -c;
prob.buc = -c;
prob.blx = -inf(size(prob.a,2),1);
prob.bux = inf(size(prob.a,2),1);
top = K.f+K.l;
prob.blx(1+K.f:K.f+K.l) = 0;

if K.q(1)>0
    for i = 1:length(K.q)
        prob.cones{i}.type = 'MSK_CT_QUAD';
        prob.cones{i}.sub  = top+1:top+K.q(i);
        top = top + K.q(i);
    end
end

if K.s(1)>0
   prob = appendSDPdata(F_struc,K,prob);
end

[r,res,solvertime] = doCall(prob,param,options);

try
    x = res.sol.itr.y;
catch
    x = nan(length(model.c),1);    
end

if options.saveduals & ~isempty(x)
    D_struc = [res.sol.itr.xx];    
    top = 1;
    for i = 1:length(K.s)
        X = zeros(K.s(i));
        n = K.s(i);
        I = find(tril(ones(n)));
        X(I) = res.sol.itr.barx(top:((top+n*(n+1)/2)-1));
        X = X + tril(X,-1)';
        D_struc = [D_struc;X(:)];
        top = top + n*(n+1)/2;
    end
else
    D_struc = [];
end

problem = MosekYALMIPError(res);

function [res,sol,solvertime] = doCall(prob,param,options)

showprogress('Calling Mosek',options.showprogress);
solvertime = clock;
if options.verbose == 0
    [res,sol] = mosekopt('minimize echo(0)',prob,param);    
else
    [res,sol] = mosekopt('minimize info',prob,param);
end
solvertime = etime(clock,solvertime);

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
        problem = 9;
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

tops = 1;
for j = 1:length(K.s)
    n = K.s(j);
    tops = [tops tops(end)+n^2];
end
top = 1+K.f+K.l+sum(K.q);
[ii,jj,kk] = find(A(top:top + sum(K.s.^2)-1,:));
cols = zeros(length(ii),1);
rows = zeros(length(ii),1);
allcol = [];
allrow = [];
allcon = [];
allvar = [];
allval = [];
for j = 1:length(K.s)    
    ind = find(ii>=tops(j) & ii<=tops(j+1)-1);
    iilocal = ii(ind)-tops(j)+1;
    col = ceil(iilocal/K.s(j));
    row = iilocal - (col-1)*K.s(j);
    allcol = [allcol col(:)'];
    allrow = [allrow row(:)'];
    allvar = [allvar jj(ind(:))'];
    allval = [allval kk(ind(:))'];
    allcon = [allcon repmat(j,1,length(col))];
end
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