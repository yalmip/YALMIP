function p = cleanuppermodel(p);

% We might have created a bilinear model from an original polynomial model.
% We should use the original model when we solve the upper bound problem.
p_bilinear = p;
p = p.originalModel;

% Remove cuts
p.F_struc(p.K.f+p.KCut.l,:)=[];
p.K.l = p.K.l - length(p.KCut.l);
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

% Remove SDP cuts
if length(p.KCut.s)>0
    starts = p.K.f+p.K.l + [1 1+cumsum((p.K.s).^2)];
    remove_these = [];
    for i = 1:length(p.KCut.s)
        j = p.KCut.s(i);
        remove_these = [remove_these;(starts(j):starts(j+1)-1)'];
    end
    p.F_struc(remove_these,:)=[];
    p.K.s(p.KCut.s) = [];
end
p.lb = p_bilinear.lb(1:length(p.c));
p.ub = p_bilinear.ub(1:length(p.c));
p.bilinears = [];