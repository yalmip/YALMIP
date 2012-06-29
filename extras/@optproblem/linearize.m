function  P = linearize(P)
%LINEARIZE Linearize constraints and objective around current solution
%
%   P = LINEARIZE(P)

P.Objective = linearize(P.Objective);
P.Constraints = linearize(P.Constraints);

