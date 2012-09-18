function [E,P] = envelope(C,x)

[aux1,aux2,aux3,p] = export(C,[],sdpsettings('solver','bmibnb'));

% Copied from bmibnb
p.high_monom_model=[];
p.originalModel = p;
[p,changed] = convert_polynomial_to_quadratic(p);
for i = 1:2
    p = presolve_bounds_from_domains(p);
    p = presolve_bounds_from_modelbounds(p);
    p = presolve_bounds_from_quadratics(p);
    p = update_eval_bounds(p);
    p = update_monomial_bounds(p);
    p = presolve_bounds_from_equalities(p);
    p = update_eval_bounds(p);
    p = update_monomial_bounds(p);
    p = update_eval_bounds(p);
    p = update_monomial_bounds(p);
end
p = compile_nonlinear_table(p);

% Copied from solvelower
p_cut = p;
if ~isempty(p.bilinears)
    p_cut.F_struc(1:p.K.f,:)=[];
    p_cut = addBilinearVariableCuts(p_cut);
    p_cut.F_struc = [p.F_struc(1:p.K.f,:);p_cut.F_struc];
end
if ~isempty(p.evalMap)
    p_cut = addEvalVariableCuts(p_cut);
end
p_cut = addMonomialCuts(p_cut);

p_cut = mergeBoundsToModel(p_cut);
if nargin > 1
    % Now project onto the variables of interest
    for i = 1:length(x)
        xi(i) = find(getvariables(x(i)) == p.used_variables);
    end
    A = -p_cut.F_struc(:,2:end);
    b = p_cut.F_struc(:,1);
    
    Akeep = A(:,xi);
    A(:,xi)=[];
    
    A = [Akeep A];
    
    Ae = A(1:p_cut.K.f,:);
    be = b(1:p_cut.K.f,:);
    A = A(1+p_cut.K.f:end,:);
    b = b(1+p_cut.K.f:end,:);
    P = Polyhedron('A',A,'b',b,'Ae',Ae,'be',be);
    P = projection(P,1:length(xi));
    
   % P = projection(polytope([Akeep A],b),1:length(xi));
    E = ismember(x,P);
else
    E = p_cut.F_struc*[1;recover(p_cut.used_variables)]>=0;
end

function p = mergeBoundsToModel(p);

A = [];
b = [];
if ~isempty(p.lb)
    A = [eye(length(p.c))];
    b = p.ub;
end
if ~isempty(p.ub)
    A = [A;-eye(length(p.c))];
    b = [b;-p.lb];
end
infbounds = find(isinf(b));
A(infbounds,:)=[];
b(infbounds)=[];
if length(b)>0
    p.F_struc = [p.F_struc(1:p.K.f,:);[b -A];p.F_struc(p.K.f+1:end,:)];
    p.K.l = p.K.l + length(b);
end







