function isherm=isessentiallyhermitian(X)
%ISHERMITIAN Check if variable is essentially Hermitian

[n,m] = size(X);
issym = 0;

if (n==m)
  isherm = norm(X-X',1)<1e-10; % Should be scaled with size maybe
end

  
  
