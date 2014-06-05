function issym=isessentiallysymmetric(X)
%ISESSENTIALLYSYMMETRIC Check if variable is essentially symmetric

% Weird name to account for MATLAB2014 and later where issymmetric exist

[n,m] = size(X);
issym = 0;

if (n==m)
  issym = norm(X-X.',1)<1e-10; % Should be scaled with size maybe
end

  
  
