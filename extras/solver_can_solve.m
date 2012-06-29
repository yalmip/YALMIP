function s = solver_can_solve(solver,problem)

% This file is currently very rudimentary, and will be much more general in
% later releases. At the moment, it is only used very limited in bmibnb, t
% osee if the upper bound solver can solve a polynomial problem without
% bilinearizing a problem

polynomial = any(problem.variabletype == 3);
sigmonial = any(problem.variabletype == 4);

s = 1;
s = s & (solver.constraint.equalities.polynomial >= polynomial);
s = s & (solver.constraint.inequalities.elementwise.polynomial >= polynomial);
s = s & (solver.objective.polynomial >= polynomial);

s = s & (solver.constraint.equalities.sigmonial >= sigmonial);
s = s & (solver.constraint.inequalities.elementwise.sigmonial >= sigmonial);
s = s & (solver.objective.sigmonial >= sigmonial);

