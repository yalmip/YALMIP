function X=flipud(X)
%FLIPUD (overloaded)

% Author Johan Löfberg 
% $Id: flipud.m,v 1.1 2006-08-10 18:00:20 joloef Exp $   

n  = X.dim(1);
m  = X.dim(2);
nm = n*m;
for i = 1:length(X.lmi_variables)+1
  X.basis(:,i) = reshape(flipud(reshape(X.basis(:,i),n,m)),nm,1);
end
% Reset info about conic terms
X.conicinfo = [0 0];

