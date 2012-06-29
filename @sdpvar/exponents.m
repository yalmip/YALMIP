function [ML,dummy] = exponents(poly,x)
%EXPONENTS Internal function to extract powers of nonlinear expression

% Author Johan Löfberg
% $Id: exponents.m,v 1.3 2006-08-11 11:48:15 joloef Exp $

mt = yalmip('monomtable');
x_lin = getvariables(poly);

x_var = getvariables(x);

ML = mt(x_lin,x_var);
if any(full(poly.basis(:,1))) %any(ismember(1,poly))
    ML = [zeros(1,length(x));ML];
end
dummy = [];