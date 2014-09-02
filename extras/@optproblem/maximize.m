function diagnostics = maximize(OptimizationProblem,Options)
%MINIMIZSE  Maximize the objective in an optimization problem
%
%   DIAGNOSTICS = MAXIMIZE(P)
%
%   Example
%
%    The following code creates an optimization problem, and then maximizes
%    the objective function 
%
%    x = sdpvar(1);P = optproblem(1 >= x >= 0, x);maximize(P)

if nargin < 2
    diagnostics = solvesdp(OptimizationProblem.Constraints,-OptimizationProblem.Objective,OptimizationProblem.Options);
else
    diagnostics = solvesdp(OptimizationProblem.Constraints,-OptimizationProblem.Objective,Options);
end
