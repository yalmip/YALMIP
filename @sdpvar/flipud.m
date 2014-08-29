function X=flipud(X)
%FLIPUD (overloaded)

n  = X.dim(1);
m  = X.dim(2);
nm = n*m;
for i = 1:length(X.lmi_variables)+1
  X.basis(:,i) = reshape(flipud(reshape(X.basis(:,i),n,m)),nm,1);
end
% Reset info about conic terms
X.conicinfo = [0 0];

