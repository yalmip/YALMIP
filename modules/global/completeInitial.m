function p = completeInitial(p)

if p.options.warmstart && any(isnan(p.x0))
  
    p_reduced = p;
    
    % Fix the given variables
    p_reduced.lb(~isnan(p.x0)) = p.x0(~isnan(p.x0));
    p_reduced.ub(~isnan(p.x0)) = p.x0(~isnan(p.x0));    
    
    % Propagate these values a bit
    p_reduced = propagate_bounds_from_equalities(p_reduced); 
    p_reduced = propagate_bounds_from_evaluations(p_reduced);
    p_reduced = propagate_bounds_from_equalities(p_reduced); 
    
    % Keep only equalities and elementwise inerqualities (more?)
    p_reduced.F_struc = p.F_struc(1:p.K.f + p.K.l,:);   
    p_reduced.K.q = 0;
    p_reduced.K.r = 0;
    p_reduced.K.p = 0;
    p_reduced.K.s = 0;
    
    % No objective
    p_reduced.c = p_reduced.c*0;
    p_reduced.Q = p_reduced.Q*0;
    
    % Kill nonlinearity information
    p_reduced.variabletype = p_reduced.variabletype*0;
    
    % Solve LP
    output = feval(p.solver.lpcall,p_reduced);
    if output.problem == 0
        p.x0 = output.Primal;
    end
end