function P = robustify(P)
%ROBUSTIFY  Derives robust counterpart.

[P.Constraints,P.Objective] = robustify(P.Constraints,P.Objective,P.Options);
