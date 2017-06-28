function p = completeInitial(p)

if p.options.usex0 && any(isnan(p.x0))
    eqKeep = ones(1,p.K.f);
    inKeep = ones(1,p.K.l);
    for i = 1:p.K.f
        if ~any(p.variabletype(find(p.F_struc(i,2:end))))
            eqKeep(i) = 1;
        else
            eqKeep(i) = 0;
        end
    end
    for i = 1:p.K.l
        if ~any(p.variabletype(find(p.F_struc(p.K.f + i,2:end))))
            inKeep(i) = 1;
        else
            inKeep(i) = 0;
        end
    end
    p_reduced = p;
    p_reduced.F_struc = [p.F_struc(find(eqKeep),:);p.F_struc(p.K.f + find(inKeep),:)];
    p_reduced.K.f = nnz(eqKeep);
    p_reduced.K.l = nnz(inKeep);
    p_reduced.K.q = 0;
    p_reduced.K.r = 0;
    p_reduced.K.p = 0;
    p_reduced.K.s = 0;
    p_reduced.c = p_reduced.c*0;
    p_reduced.Q = p_reduced.Q*0;
    p_reduced.variabletype = p_reduced.variabletype*0;
    p_reduced.lb(~isnan(p.x0)) = p.x0(~isnan(p.x0));
    p_reduced.ub(~isnan(p.x0)) = p.x0(~isnan(p.x0));
    output = feval(p.solver.lpcall,p_reduced);
    if output.problem == 0
        p.x0 = output.Primal;
    end
end