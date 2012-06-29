function Z=clean(Z,tol)
%CLEAN Remove terms with small coefficients
%
% Z = clean(X,tol) removes all variables with a coefficient smaller than tol
%
%

% Author Johan Löfberg
% $Id: clean.m,v 1.1 2005-12-19 15:27:31 joloef Exp $

if nargin == 1
    tol = 0;
end
basis_real = real(Z);
basis_imag = imag(Z);
basis_real(abs(basis_real)<tol) = 0;
basis_imag(abs(basis_imag)<tol) = 0;
Z = basis_real + sqrt(-1)*basis_imag;
