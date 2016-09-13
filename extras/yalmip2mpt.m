function Matrices = yalmip2mpt(interfacedata)

F_struc = interfacedata.F_struc;
c       = interfacedata.c;
f       = interfacedata.f;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;

% We have three groups of variables.
% parametric, free, and binary free
param_var  = interfacedata.parametric_variables;
free_var   = setdiff(1:length(c),param_var);
binary_var = interfacedata.binary_variables;

% These rows of the solution should be returned
if isempty(interfacedata.requested_variables)
    Matrices.requested_variables = 1:length(free_var);
else
    Matrices.requested_variables = [];
    for i = 1:length(interfacedata.requested_variables)
        Matrices.requested_variables = [Matrices.requested_variables  find(ismember(free_var,interfacedata.requested_variables(i)))];
    end
end

% Extract data for free and parametric variables
Matrices.F  = (2*Q(param_var,free_var));
Matrices.Y  = (Q(param_var,param_var));
Matrices.H  = full((2*Q(free_var,free_var)));
Matrices.G  = (-F_struc(1+K.f:end,1+free_var));
Matrices.E  = (F_struc(1+K.f:end,1+param_var));
Matrices.W  = (F_struc(1+K.f:end,1));
Matrices.Cf = (c(free_var)');
Matrices.Cx = (c(param_var)');
Matrices.Cc = (f);

% Equality constraints (YALMIP specific)
Matrices.Aeq = -F_struc(1:K.f,1+free_var);
Matrices.Beq = -F_struc(1:K.f,1+param_var);
Matrices.beq = F_struc(1:K.f,1);

% This is not dealth with in YALMIP
Matrices.bndA=[];
Matrices.bndb=[];

% Save this silly info also
Matrices.nu = length(free_var);
Matrices.nx = length(param_var);

if nnz(Matrices.H)==0 & nnz(Matrices.Y)==0 && ~isequal(interfacedata.solver.tag,'POP')
    % Whoops, it's an LP and MPT uses the name H for linear cost also...
    Matrices.H = c(free_var)';    
    if nnz(Matrices.F) > 0
        % Latest MPT does not fullify D
        Matrices.D = full(Matrices.F');
    end
    Matrices.F = Matrices.Cx;
    Matrices.qp = 0;
else
    Matrices.qp = 1;
end

if ~isempty(lb)
    Matrices.lb = [lb(free_var);lb(param_var)];
else
    Matrices.lb = repmat(-inf,nx+nu,1);
end
if ~isempty(ub)
    Matrices.ub = [ub(free_var);ub(param_var)];
else
    Matrices.ub = repmat(inf,nx+nu,1);
end 

Matrices.param_var  = interfacedata.parametric_variables;
Matrices.free_var   = setdiff(1:length(c),param_var);
Matrices.binary_var_index = find(ismember(free_var,binary_var));

% A transformation to go to original variables when 
% equality constraints have been removed
Matrices.getback = [];