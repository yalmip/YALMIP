function [yesno,diagnostic] = isfeasible(OptimizationProblem,Options)
%FEASIBLE  Check feasibility by solving feasibility problem
%
%   YESNO = isfeasible(P[,Options])
%
%   Example
%
%    The following code creates an optimization problem, and then checks
%    feasibility
%
%    x = sdpvar(1);P = optproblem(x >= 0, x^2);isfeasible(P)

OptimizationProblem.Options.verbose = 0;
if nargin < 2
    diagnostics = solvesdp(OptimizationProblem.Constraints,[],OptimizationProblem.Options);
    yesno = diagnostics.problem ~=1;
else
    Options.verbose = 0;
    diagnostics = solvesdp(OptimizationProblem.Constraints,[],OptimizationProblem.Objective,Options);
    yesno = diagnostics.problem ~=1;
end
