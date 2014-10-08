function [p,c,v] = matrixpolynomial(x,n,dmax,dmin)
%MATRIXPOLYNOMIAL Creates parameterized polynomial
%
% [p,c,v] = matrixpolynomial(x,n,dmax,dmin)
%
% MATRIXPOLYNOMIAL is a quick way to define a parameterized polynomial 
% p=sum C_i v_i(x) with all monomials of dmin <= degree(p,x) <= dmax. The
% coefficients in the polynomial are C_i while v is the monomial basis.
%
% Example:
%
% Paramterized quartic 2x2 matrix
%  x = sdpvar(2,1);
%  p = matrixpolynomial(x,2,4);
%
% See also MONOLIST, COEFFICIENTS

if (length(dmax) > 1) && (length(dmax) ~= length(x))
    error('Dimension mismatch: The third argument should be the max degree for each variable, or a sclar');
end

if nargin > 3
    if (length(dmin) > 1) && (length(dmin) ~= length(x))
        error('Dimension mismatch: The third argument should be the max degree for each variable, or a sclar');
    end
end

if any(dmax < 0)
    error('Only non-negative polynomial degrees possible')
end

if nargin<4
    dmin = 0;
end

if any(dmin > dmax)
    error('Fourth argument (dmin) should not be larger than second argument (dmax)');
end

if any(dmin < 0)
    error('Only non-negative polynomial degrees possible')
end

if length(n)==1
    n = [n n];
end

v = monolist(x,dmax);

if dmin <= dmax & dmin>0
    s = nchoosek(length(x) + dmin-1,dmin-1);
    v = extsubsref(v,s+1:length(v));
end

p = 0;
c = [];
for i = 1:length(v)
    Ci = sdpvar(n(1),n(2));
    c = [c reshape(Ci,[],1)];
end
p = reshape(c*v,n(1),n(2));

