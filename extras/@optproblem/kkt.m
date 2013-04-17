function [KKTConstraints, details] = kkt(P,x,ops)
%KKT Create KKT system for optimization system P with parametric variables x
%
% [KKTConstraints, details] = kkt(Constraints,Objective,parameters,options)

if nargin < 2
    [KKTConstraints, details] = kkt(P.Constraints,P.Objective);
elseif nargin < 3
    [KKTConstraints, details] = kkt(P.Constraints,P.Objective,x);
else
    [KKTConstraints, details] = kkt(P.Constraints,P.Objective,x,ops);
end
