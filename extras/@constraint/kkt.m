function [KKTConstraints, details] = kkt(F,h,parametricVariables,ops);
%KKT Create KKT system
%
% [KKTConstraints, details] = kkt(Constraints,Objective,parameters)

% Author Johan Löfberg
% $Id: kkt.m,v 1.2 2010-02-08 13:06:11 joloef Exp $

if nargin == 2
    [KKTConstraints, details] = kkt(set(F),h);
elseif nargin==3
    [KKTConstraints, details] = kkt(set(F),h,parametricVariables);    
else
    [KKTConstraints, details] = kkt(set(F),h,parametricVariables,ops);    
end

