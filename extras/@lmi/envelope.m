function [E,P] = envelope(C,x)

[aux1,aux2,aux3,p] = export(C,[],sdpsettings('solver','bmibnb'));

% Copied from bmibnb
p.high_monom_model=[];
p0 = p;
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
if any(p0.variabletype==3)
    monomials = find(p0.variabletype==3);
    for i = 1:length(monomials)
        monom_index = monomials(i);
        monom = p0.monomtable(monomials(i),:);
        monom_variable = find(monom);
        if length(monom_variable)==1
            n = monom(monom_variable);
            if ~even(n)
                if p.lb(monom_variable)<0 &  p.ub(monom_variable)>0
                    L = p.lb(monom_variable);
                    U = p.ub(monom_variable);
                    
                    % Tangent at x = lower bound
                    if L<=0
                        p_cut.F_struc(end+1,1) = L^n-n*L^n;
                        p_cut.F_struc(end,1+monom_index)=-1;
                        p_cut.F_struc(end,1+monom_variable)=n*L^(n-1);
                        p_cut.K.l = p_cut.K.l+1;
                    end
                     
                    % Tangent at x = upper bound
                    if U >= 0
                        p_cut.F_struc(end+1,1) = -(U^n-n*U^n);
                        p_cut.F_struc(end,1+monom_index)= 1;
                        p_cut.F_struc(end,1+monom_variable)=-n*U^(n-1);
                        p_cut.K.l = p_cut.K.l+1;
                    end
                    
                    % Line between lower bound and tangent intersection
                    r = zeros(1,n+1);r(1)=n-1;r(2)=-L*n;r(end)=L^n;
                    r = roots(r);
                    r = r(min(find(r==real(r))));
                    if r >= U
                        r = U;
                        fprim = (U^n-L^n)/(U-L);
                    else
                        fprim = n*r^(n-1);
                    end
                    p_cut.F_struc(end+1,1) = -L^n+L*fprim;
                    p_cut.F_struc(end,1+monom_index)=1;
                    p_cut.F_struc(end,1+monom_variable)=-fprim;
                    p_cut.K.l = p_cut.K.l+1;
                    
                    % Line between upper bound and tangent intersection
                    r = zeros(1,n+1);r(1)=n-1;r(2)=-U*n;r(end)=U^n;
                    r = roots(r);
                    r = r(min(find(r==real(r))));
                    if r <= L
                        r = L;
                        fprim = (U^n-L^n)/(U-L);
                    else
                        fprim = n*r^(n-1);
                    end
                    p_cut.F_struc(end+1,1) = U^n-U*fprim;
                    p_cut.F_struc(end,1+monom_index)=-1;
                    p_cut.F_struc(end,1+monom_variable)=fprim;
                    p_cut.K.l = p_cut.K.l+1;
                end
            end
        end
    end
end

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
    P =  projection(P,1:length(xi));
    
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







