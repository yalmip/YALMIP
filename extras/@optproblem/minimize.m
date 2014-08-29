function diagnostics = minimize(OptimizationProblem,Options)
%MINIMIZSE  Minimize the objective in an optimization problem
%
%   DIAGNOSTICS = MINIMIZE(P,Options)
%
%   Example
%
%    The following code creates an optimization problem, and then minimizes
%    the objective function 
%
%    x = sdpvar(1);P = optproblem(x >= 0, x^2);minimize(P)

if nargin < 2
    diagnostics = solvesdp(OptimizationProblem.Constraints,OptimizationProblem.Objective,OptimizationProblem.Options);
else
    diagnostics = solvesdp(OptimizationProblem.Constraints,OptimizationProblem.Objective,Options);
end
