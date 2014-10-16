function [model,nonlinearremain] = yalmip2cplex(interfacedata)

% Author Johan Löfberg

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;

ub      = interfacedata.ub;
lb      = interfacedata.lb;
H       = 2*interfacedata.Q;
if ~isempty(semicont_variables)
    [F_struc,K,lb,ub,semicont_variables] = extractSemiContBounds(F_struc,K,lb,ub,semicont_variables);
end
[F_struc,K,c,H,ub,lb,x0,Qi,Li,ri] = append_normalized_socp(F_struc,K,c,H,ub,lb,x0);

nonlinearremain=0;
if any(interfacedata.variabletype)
    nonlinearremain=1;    
end

% Notation used
f = c;
if ~isempty(F_struc)
    Aineq = -F_struc(:,2:end);
    bineq = F_struc(:,1);
else
    Aineq = [];
    bineq = [];
end

% CPLEX assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
if ~isempty(semicont_variables)
    [NegativeSemiVar,H,f,Aineq,lb,ub,semicont_variables] = negateNegativeSemiContVariables(H,f,Aineq,lb,ub,semicont_variables,[]);
else
    NegativeSemiVar = [];
end

if K.f > 0
    Aeq = Aineq(1:K.f,:);
    beq = bineq(1:K.f);
    % Code around performance bugs in older version of MATLAB
    [ii,jj,ss] = find(Aineq);keeps = ii>K.f;
    Aineq = sparse(ii(keeps)-K.f,jj(keeps),ss(keeps),size(Aineq,1)-K.f,size(Aineq,2));
    [ii,jj,ss] = find(bineq);keeps = ii>K.f;
    bineq = sparse(ii(keeps)-K.f,jj(keeps),ss(keeps),length(bineq)-K.f,1);
else
    Aeq = [];
    beq = [];
end

if all(isinf(lb))
    lb = [];
end
if all(isinf(ub))
    ub = [];
end

ctype = char(ones(length(f),1)*67);
if ~isempty(integer_variables)
    ctype(setdiff(integer_variables,semicont_variables)) = 'I';
end
ctype(binary_variables)  = 'B';  % Should not happen except from bmibnb
if ~isempty(semicont_variables)
    ctype(setdiff(semicont_variables,integer_variables)) = 'S';  % Should not happen except from bmibnb
    ctype(intersect(semicont_variables,integer_variables)) = 'N';
end

options.cplex.simplex.display = options.verbose;
options.cplex.mip.display = options.verbose;
options.cplex.barrier.display = options.verbose;
if str2num(interfacedata.solver.subversion)>=12.3
    if options.verbose
        options.cplex.diagnostics = 'on';
    else
        options.cplex.diagnostics = 'off';
    end
end

if ~isempty(K.sos.type)
    for i = 1:length(K.sos.weight)
        K.sos.weight{i} = full(K.sos.weight{i}(:));
    end
    if length(K.sos.weight)==1
        K.sos.weight = K.sos.weight{1};
        K.sos.variables = K.sos.variables{1};
    end
end

% Shift the objective (constant term can not be set in CPLEX?)
if isfield(options.cplex,'mip.tolerances.lowercutoff')
    options.cplex.mip.tolerances.lowercutoff = options.cplex.mip.tolerances.lowercutoff-interfacedata.f;
    options.cplex.mip.tolerances.uppercutoff = options.cplex.mip.tolerances.uppercutoff-interfacedata.f;
end

if nnz(H)==0
    H = [];
end

% Bug in cplex 12.1
if length(K.q) == 1
    if isa(Qi,'cell')    
        Qi = Qi{1};
    end
end 

model.H = H;
model.f = f;
model.Aineq = Aineq;
model.bineq = bineq;
model.Aeq = Aeq;
model.beq = beq;
model.lb = lb;
model.ub = ub;
model.x0 = x0;
model.options = options.cplex;
model.verbose = options.verbose;
model.integer_variables = integer_variables;
model.binary_variables = binary_variables;
model.semicont_variables = semicont_variables;
model.K = K;
model.ctype = ctype;
model.Qi=Qi;
model.Li=Li;
model.ri=ri;
model.NegativeSemiVar=NegativeSemiVar;

