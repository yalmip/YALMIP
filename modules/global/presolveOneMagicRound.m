function p = presolveOneMagicRound(p);

p = presolve_bounds_from_domains(p);
p = presolve_bounds_from_modelbounds(p);
p = propagate_bounds_from_convex_quadratic_ball(p);
p = propagate_bounds_from_evaluations(p);
p = update_monomial_bounds(p);
p = propagate_bounds_from_equalities(p);
p = presolve_bounds_from_inequalities(p);
p = propagate_bounds_from_evaluations(p);
p = update_monomial_bounds(p);
p = propagate_bounds_from_evaluations(p);
p = update_monomial_bounds(p);
p = propagate_bounds_from_arbitrary_quadratics(p);