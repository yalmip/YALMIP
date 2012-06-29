function res = isinterval(Y)
%ISINTERVAL (overloaded)

% Author Johan Löfberg 
% $Id: isinterval.m,v 1.1 2006-12-14 13:20:36 joloef Exp $   

res = isa(Y.basis,'intval');
