function issym=issymmetric(X)
%ISSYMMETRIC Check if variable is symmetric

% Author Johan Löfberg 
% $Id: issymmetric.m,v 1.3 2008-04-02 19:06:54 joloef Exp $

[n,m] = size(X);
issym = 0;

if (n==m)
  issym = norm(X-X.',1)<1e-10; % Should be scaled with size maybe
end

  
  
