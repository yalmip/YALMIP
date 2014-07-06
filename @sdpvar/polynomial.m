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

if any(dmax < 0)
    error('Only non-negative polynomial degrees possible')
end

if nargin<3
    dmin = 0;
end

if any(dmin > dmax)
    error('Third argument (dmin) should not be larger than second argument (dmax)');
end

if any(dmin < 0)
    error('Only non-negative polynomial degrees possible')
end

v = monolist(x,dmax);

if dmin <= dmax & dmin>0
    s = nchoosek(length(x) + dmin-1,dmin-1);
    v = extsubsref(v,s+1:length(v));
end
c = sdpvar(length(v),1);
p = c'*v;