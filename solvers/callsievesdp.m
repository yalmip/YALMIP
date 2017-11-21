function [Constraints,Objective] = callsievesdp(Constraints,Objective,options)

% Export in mosek format and reduce model
[model,recovery,diagnostics,internal] = export(Constraints, Objective, sdpsettings('solver','mosek'));
y_original = recover(recovery.used_variables);
probOriginal = model.prob;
[probReduced,info] = SieveSDP(probOriginal);

% Now go back to YALMIP model and solve the problem with the solver
% specified (quick hack, will be done in purely numerical format in real
% release)
[A, b, c, K] = convert_mosek2sedumi(probReduced);
[Constraint_r, Objective_r, y_reduced] = sedumi2yalmip(A,b,c,K);

if nargin < 3
    options = [];
end
optimize(Constraint_r,Objective_r,options);

% Problem: How do we go back from the variables in the reduced model to the
% variables in the original 
%assign(y_original, something value(y_reduced))
