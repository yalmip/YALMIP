function diagnostics = feasible(OptimizationProblem,Options)
%FEASIBLE  Computes feasible solution of in an optimization problem
%
%   DIAGNOSTICS = FEASIBLE(P,Options)
%
%   Example
%
%    The following code creates an optimization problem, and then finds a
%    feasible solution
%
%    x = sdpvar(1);P = optproblem(1 >= x >= 0, x^2);feasible(P)

if nargin < 2
    diagnostics = solvesdp(OptimizationProblem.Constraints,[],OptimizationProblem.Options);
else
    diagnostics = solvesdp(OptimizationProblem.Constraints,[],Options);
end
