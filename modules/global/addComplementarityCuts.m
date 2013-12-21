function p = addComplementarityCuts(p)

if size(p.complementary,1)>0
    p_cut = emptyNumericalModel;
    for i = 1:size(p.complementary,1)
        variables = p.complementary(i,:);
        if p.lb(variables(1))>0 & p.lb(variables(2))>0
            p_cut.F_struc(end+1,1)=p.ub(variables(2));
            p_cut.F_struc(end,1+variables(1))=-p.ub(variables(2))/p.ub(variables(1));
            p_cut.F_struc(end,1+variables(2))=-1;
            p_cut.K.l = p_cut.K.l + 1;
        end
    end
    p = mergeNumericalModels(p,p_cut);
end