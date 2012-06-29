function [p,c,v] = polynomial(x,dmax,dmin)
%POLYNOMIAL Creates parameterized polynomial
%
% [p,c,v] = polynomial(x,dmax,dmin)
%
% POLYNOMIAL is a quick way to define a parameterized polynomial p=c'*v,
% with all monomials of dmin <= degree(p,x) <= dmax. The coefficients in
% the polynomial are c while v is the monomial basis.
%
% Example:
%
% Paramterized quartic
%  x = sdpvar(2,1);
%  p = polynomial(x,4);
%
% See also MONOLIST, COEFFICIENTS

% Author Johan Löfberg 
% $Id: polynomial.m,v 1.3 2006-08-16 15:57:36 joloef Exp $

if nargin<3
    dmin = 0;
end

v = monolist(x,dmax);
if dmin <= dmax & dmin>0
    s = nchoosek(length(x) + dmin-1,dmin-1);
    v = extsubsref(v,s+1:length(v));
end
c = sdpvar(length(v),1);
p = c'*v;