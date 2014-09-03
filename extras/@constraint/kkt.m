function [KKTConstraints, details] = kkt(F,h,parametricVariables,ops);
%KKT Create KKT system
%
% [KKTConstraints, details] = kkt(Constraints,Objective,parameters)

if nargin == 2
    [KKTConstraints, details] = kkt(lmi(F),h);
elseif nargin==3
    [KKTConstraints, details] = kkt(lmi(F),h,parametricVariables);    
else
    [KKTConstraints, details] = kkt(lmi(F),h,parametricVariables,ops);    
end

