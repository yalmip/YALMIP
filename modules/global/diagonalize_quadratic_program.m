function [pnew,timing,x_min,upper] = diagonalize_quadratic_program(p,timing,x_min,upper)
pnew = p;
pnew.V = [];
pnew.diagonalized = 0;

if ~p.options.bmibnb.diagonalize
    return
end

% Preprocessing has moved nonliner constraints
% to bounds, so don't be fooled...
if p.originallyNonlinearConstraints
    return
end

% No quadratic terms
if all(p.variabletype == 0)
    return
end

% Any polynomial terms or simple linear by some reason
if any(p.variabletype > 2) & ~all(p.variabletype == 0) | p.options.bmibnb.diagonalize==0
    return
end

if ~isempty(p.evalVariables)
    return
end

if ~isempty(p.binary_variables) | ~isempty(p.integer_variables)
    return
end

if isempty(p.F_struc)
    % Bound-constrained QP is faster in original
    % form as we can exploit concavity
    return
end

nonlinear = find(p.variabletype > 0);
linear = find(p.variabletype == 0);
if ~isempty(p.F_struc)
    % Nonlinear terms in constraints
    if nnz(p.F_struc(:,1 + nonlinear))> 0
        return
    end
end

% Find quadratic and linear terms
used_in_c = find(p.c);
quadraticterms = used_in_c(find(ismember(used_in_c,nonlinear)));
Q = zeros(length(p.c),length(p.c));
if ~isempty(quadraticterms)
    usedinquadratic = zeros(1,length(p.c));
    for i = 1:length(quadraticterms)
        Qij = p.c(quadraticterms(i));
        power_index = find(p.monomtable(quadraticterms(i),:));
        if length(power_index) == 1
            Q(power_index,power_index) = Qij;
        else
            Q(power_index(1),power_index(2)) = Qij/2;
            Q(power_index(2),power_index(1)) = Qij/2;
        end       
    end
end
Qlin = Q(linear,linear);
clin = p.c(linear);

if nnz(Qlin-diag(diag(Qlin)))==0
    % already diagonal
    return
end

% Decompose Q
[V,D] = eig(full(Qlin));
V = real(V);
D = real(D);
V(abs(V)<1e-11) = 0;
D(abs(D)<1e-11) = 0;
lb = p.lb(linear);
ub = p.ub(linear);
Z = V';
newub = sum([Z>0].*Z.*repmat(ub',length(Z),1),2)+sum([Z<0].*Z.*repmat(lb',length(Z),1),2);
newlb = sum([Z>0].*Z.*repmat(lb',length(Z),1),2)+sum([Z<0].*Z.*repmat(ub',length(Z),1),2);
newub(isnan(newub)) = inf;
newlb(isnan(newlb)) = -inf;

% Create new problem
clin = V'*clin;

n = length(linear);
pnew.original_linear = linear;
pnew.original_n = length(p.c);
pnew.V = V;
pnew.c = [clin;diag(D)];
pnew.Q = spalloc(2*n,2*n,0);

% find constraint polytope
if size(p.F_struc,1)>0
    A = -p.F_struc(:,1 + linear);
    b = p.F_struc(:,1);
    pnew.F_struc = [b -A*V zeros(length(b),n)];
    Abounds = [V;-V];
    bbounds = [ub;-lb];
    keep = find(~isinf(bbounds));
    if ~isempty(keep)
        pnew.F_struc = [pnew.F_struc;bbounds(keep) -Abounds(keep,:) zeros(length(keep),n)];       
        pnew.K.l = pnew.K.l + length(keep);
    end
end

pnew.variabletype = [zeros(1,n) ones(1,n)*2];
pnew.monomtable = [eye(n);2*eye(n)];
pnew.monomtable(2*n,2*n) = 0;
pnew.lb =-inf(2*n,1);
pnew.ub = inf(2*n,1);
pnew.lb(1:n) = newlb;
pnew.ub(1:n) = newub;
pnew.lb(n+1:end) = 0;
pnew.ub(n+1:end) = max(newlb.^2,newub.^2);
if length(p.x0)>0
    pnew.x0 = [(V'*p.x0(p.linears));(V'*p.x0(p.linears)).^2];
end
x_min = [V'*x_min(p.linears);(V'*x_min(p.linears)).^2];
pnew.diagonalized = 1;

% Ok ,we now got a new diagonalized model
% Sometimes it performs much worse though
% Make trial calls to judge it
if p.options.bmibnb.diagonalize == -1
    tstart = tic; 
    % If the model lacks bounds, these relaxations will
    % fail at this point, as bound propagation hasn't been performed yet.
    % Add fake bounds and see which performs best on fake model       
    p.socpcuts.F_struc = [];
    pnew.socpcuts.F_struc = [];
    p.delayedconvex = [];
    pnew.delayedconvex = [];
    p.upper=inf;
    pnew.upper=inf;
    [output1,cost1,~,timing] = solvelower(fakeBounds(defineQuadratics(fakeLowerModel(pnew))),p.options,p.solver.lowersolver.call,[],[],timing);
    [output2,cost2,~,timing] = solvelower(fakeBounds(fakeLowerModel(p)),p.options,p.solver.lowersolver.call,[],[],timing);
    [upper_cost1,x_u1] = quickCost(pnew,output1.Primal);
    [upper_cost2,x_u2] = quickCost(p,output2.Primal);
    timing.lowersolve = timing.lowersolve + toc(tstart);
    if output1.problem && output2.problem
        % both had issues. Let's assume diagonal is better
        if p.options.verbose>=0;display('* -Diagonalized QP');end
    elseif output1.problem && ~output2.problem
        % new model has issues. use old
        pnew = p;
        pnew.V = [];
        pnew.diagonalized = 0;
        upper = upper_cost2;x_min = x_u2;
        if p.options.verbose>=0;display('* -Diagonalized QP but switched back to original.');end
    elseif ~output1.problem && output2.problem
        % old has issue, so keep new
        upper = upper_cost1;x_min = x_u1;
        if p.options.verbose>=0;display('* -Diagonalized QP.');end
    else
        % Both solved, so which had the best lower bound
        if cost1 < cost2
            % New model weaker, keep old. 
            % Take best upper though!
            x_min = x_u2;
            upper = upper_cost2;
            if upper_cost1 < upper
                upper = upper_cost1;
                x_min(p.linears)=V*x_u1(pnew.linears);
            end
            pnew = p;
            pnew.V = [];
            pnew.diagonalized = 0;
            if p.options.verbose>=0;display('* -Diagonalized QP but switched back to original.');end
        else
            % New model better
            % Take best upper though!
            x_min = x_u1;
            upper = upper_cost1;
            if upper_cost2 < upper
                upper = upper_cost2;
                x_min(pnew.linears)=V'*x_u2(p.linears);
            end
            if p.options.verbose>=0;display('* -Diagonalized QP.');end
        end
    end
else
	if p.options.verbose>=0;display('* -Diagonalized QP.');end    
end
       
function pnew = fakeLowerModel(p)
pnew = p;
pnew.complementary=[];
pnew.shiftedQP=[];
pnew.originalModel.variabletype=[];
pnew.concavityEqualities=[];
pnew.EqualityConstraintState=[];
pnew.EqualityConstraintState=[];
pnew.InequalityConstraintState=[];

function [q,z] = quickCost(p,x)
q = 0;
z=[];
for i = 1:size(p.monomtable,1)
    j = find(p.monomtable(i,:));
    v = x(j);
    w = p.monomtable(i,j);
   	z =[z;prod(v(:).^w(:))];
    q = q + p.c(i)*z(i);
end

function p = defineQuadratics(p)
n = length(p.c)/2;
p.Quadratics = n+1:2*n;
p.nonlinears = n+1:2*n;
p.QuadraticsList = [zeros(n,2);(1:n)' (1:n)'];
p.bilinears = [(n+1:2*n)' (1:n)' (1:n)'];
%p.solver.lowersolver.constraint.inequalities.secondordercone.linear=0;

function p = fakeBounds(p)
U = 1e3;
p.ub(p.linears) = min(p.ub(p.linears),U);
p.lb(p.linears) = max(p.lb(p.linears),-U);
p.ub(p.Quadratics) = min(p.ub(p.Quadratics),U^2);
p.lb(p.Quadratics) = max(p.lb(p.Quadratics),-U^2);