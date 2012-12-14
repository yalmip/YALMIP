function [p,E] = poly (A, r)
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

p{n+1} = 1;
E = [];
for k = n: -1: 1,
  M = extsubsref(A,k+1:n,k+1:n);
  R = extsubsref(A,k,k+1:n);
  S = extsubsref(A,k+1:n,k);
  ci = [1 -extsubsref(A,k,k)];
  Ri = R;
  for i = k+1:n
    ci = [ci -R*S];
    if (i < n) 
        if isa(R,'double')
            R = R*M;
        else
         Ri = R;
         R = sdpvar(size(Ri,1),size(Ri,2),'full');
         E = [E, R == Ri*M]; 
        end
    end
  end
  C = toeplitz (ci, [1 zeros(1,n-k)]);
  p{k} = sdpvar(size(C,1),1);
  if isa(p{k+1},'double')
      p{k} = C*p{k+1};
  else
       ptemp = C * p{k+1};
       for i = 1:length(ptemp)
           if isa(extsubsref(ptemp,i),'double')
               ttt(i,1)=0;
           else
           ttt(i,1) = 1;
           end
       end
       g = getbase(ptemp);
       E = [E, p{k}.*ttt == (C * p{k+1}).*ttt];
       p{k} = p{k}.*ttt + (1-ttt).*ptemp;
  end
end
p = [p{1}]';

if nargin >= 2,
    p = extsubsref(p,n+1-r);
end
