function E = envelope(C,x)

[aux1,aux2,aux3,p] = export(C,[],sdpsettings('solver','bmibnb'));

% Copied from bmibnb
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
p = compile_nonlinear_table(p);
[p,x_min,upper] = initializesolution(p);

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

% Now project onto the variables of interest
for i = 1:length(x)
    xi(i) = find(getvariables(x(i)) == p.used_variables);
end
A = -p_cut.F_struc(:,2:end);
b = p_cut.F_struc(:,1);

Akeep = A(:,xi);
A(:,xi)=[];
P = projection(polytope([Akeep A],b),1:length(xi));
E = ismember(x,P);



