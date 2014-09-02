function X=diag(X,k)
%DIAG (overloaded)

n = X.dim(1);
m = X.dim(2);

if nargin == 1
    k = 0;
end

if min([n m])==1
    Y = X;
	diagX = diag(X.basis(:,1),k);
	Y.basis = diagX(:);
	for i = 1:length(X.lmi_variables)
		diagX = diag(X.basis(:,i+1),k);
		Y.basis(:,i+1) = diagX(:);
    end
	Y.lmi_variables = X.lmi_variables;
	Y.dim(1) = size(diagX,1);
	Y.dim(2) = size(diagX,2);
    X=Y;
else
    index = diag(reshape(1:n*m,n,m),k);
    X.basis = X.basis(index,:);
	X.dim(1) = length(index);
	X.dim(2) = 1;
    X = clean(X);
end
% Reset info about conic terms
if isa(X,'sdpvar')
    X.conicinfo = [0 0];
end