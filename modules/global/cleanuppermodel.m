function p = cleanuppermodel(p)

% We might have created a bilinear model from an original polynomial model.
% We should use the original model when we solve the upper bound problem.
p_bilinear = p;
p = p.originalModel;

p = removeCuts(p);

n_start = length(p.c);

% Quadratic mode, and quadratic aware solver?
bilinear_variables = find(p.variabletype == 1 | p.variabletype == 2);
if ~isempty(bilinear_variables)
    used_in_c = find(p.c);
    quadraticterms = used_in_c(find(ismember(used_in_c,bilinear_variables)));
    if ~isempty(quadraticterms) & p.solver.uppersolver.objective.quadratic.nonconvex
        usedinquadratic = zeros(1,length(p.c));
        for i = 1:length(quadraticterms)
            Qij = p.c(quadraticterms(i));
            power_index = find(p.monomtable(quadraticterms(i),:));
            if length(power_index) == 1
                p.Q(power_index,power_index) = Qij;
            else
                p.Q(power_index(1),power_index(2)) = Qij/2;
                p.Q(power_index(2),power_index(1)) = Qij/2;
            end
            p.c(quadraticterms(i)) = 0;
        end
    end
end

p.lb = p_bilinear.lb(1:length(p.c));
p.ub = p_bilinear.ub(1:length(p.c));
p.bilinears = [];