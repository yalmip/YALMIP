function diagnostics = optimize(OptimizationProblem,Options)
%OPTIMIZE  Minimize the objective in an optimization problem
%
%   DIAGNOSTICS = OPTIMIZE(P)
%
%   Example
%
%    The following code creates an optimization problem, and then minimizes
%    the objective function
%
%    x = sdpvar(1);P = optproblem(x >= 0, x^2);optimize(P)

if nargin < 2
    diagnostics = solvesdp(OptimizationProblem.Constraints,OptimizationProblem.Objective,OptimizationProblem.Options);
else
    diagnostics = solvesdp(OptimizationProblem.Constraints,OptimizationProblem.Objective,Options);
end
