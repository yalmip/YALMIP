function p = poly (A, r)
% POLY (Overloaded) 
%
% Computes the coefficients of the characteristic polynomial, det (sI-A),
% of a square matrix, A, using Berkowitz algorithm. 
%
% p = poly(A)
%
% Stuart J. Berkowitz, "On computing determinat in small parallel time
% using a small number of processors. Information Processing Letters,
% 18(3):147-150, 1984.
%
% See also DET

% Author Anders Helmersson, Johan Löfberg
% $Id: poly.m,v 1.2 2006-10-24 12:02:04 joloef Exp $%

[n, m] = size (A);
if n ~= m, error ('A must be square'); end;

p = 1;

for k = n: -1: 1,
  M = extsubsref(A,k+1:n,k+1:n);
  R = extsubsref(A,k,k+1:n);
  S = extsubsref(A,k+1:n,k);
  ci = [1 -extsubsref(A,k,k)];
  for i = k+1:n
    ci = [ci -R*S];
    if (i < n) R = R*M; end
  end
  C = toeplitz (ci, [1 zeros(1,n-k)]);
  p = C * p;
end
p = p';

if nargin >= 2, p = p(n+1-r); end
