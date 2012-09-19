function p = presolveOneMagicRound(p);

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
