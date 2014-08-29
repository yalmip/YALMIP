function [KKTConstraints, details] = kkt(F,h,parametricVariables,ops);
%KKT Create KKT system
%
% [KKTConstraints, details] = kkt(Constraints,Objective,parameters)

if nargin == 2
    [KKTConstraints, details] = kkt(set(F),h);
elseif nargin==3
    [KKTConstraints, details] = kkt(set(F),h,parametricVariables);    
else
    [KKTConstraints, details] = kkt(set(F),h,parametricVariables,ops);    
end

