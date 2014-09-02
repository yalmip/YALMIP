function [p,c] = polynomial(x,dmax,dmin)
%POLYNOMIAL Creates parameterized polynomial
%
% [p,c] = polynomial(x,dmax,dmin)
%
% POLYNOMIAL is a quick way to define a parameterized polynomial, with all
% monomials of dmin <= degree(p,x) <= dmax.
%
% Example:
%
% Paramterized quartic
%  x = sdpvar(2,1);
%  [p,c] = polynomial(x,4);
%
% See also MONOLIST, COEFFICIENTS

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