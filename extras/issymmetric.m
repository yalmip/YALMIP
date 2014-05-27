function issym=issymmetric(X)
%ISSYMMETRIC Check if variable is symmetric

[n,m] = size(X);
issym = 0;

if (n==m)
  issym = norm(X-X.',1)<1e-10; % Should be scaled with size maybe
end

  
  
