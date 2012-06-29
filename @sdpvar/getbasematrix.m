function Q=getbasematrix(X,ind)
%GETBASEMATRIX Internal function to extract basematrix for variable IND

% Author Johan Löfberg 
% $Id: getbasematrix.m,v 1.3 2006-07-26 20:17:58 joloef Exp $  

if ind==0
  base = X.basis(:,1);
  Q = reshape(base,X.dim(1),X.dim(2));
  return;
end

here = find(X.lmi_variables==ind);
if isempty(here)
  Q = sparse(X.dim(1),X.dim(2));
else
  base = X.basis(:,here+1);
  Q = reshape(base,X.dim(1),X.dim(2));
end


