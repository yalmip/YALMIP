function X=tril(X,r)
%TRIL (overloaded)

if nargin==1
    r = 0;
end

n = X.dim(1);
m = X.dim(2);

ind = reshape(1:n*m,n,m);
ind = tril(ind,r);
ind = find(ind==0);

[i,j,k] = find(X.basis);
zero_these = find(ismember(i,ind));
i(zero_these)=[];
j(zero_these)=[];
k(zero_these)=[];
X.basis = sparse(i,j,k,size(X.basis,1),size(X.basis,2));
X.conicinfo = [0 0];
X = clean(X);
return

trilX = tril(reshape(X.basis(:,1),n,m),r);
Y.basis = trilX(:);

j = 1;
for i = 1:length(x_lmi_variables)
	trilX = tril(reshape(X.basis(:,i+1),n,m),r);
	if (norm(trilX,inf)>0)
		Y.basis(:,j+1) = trilX(:);
		lmi_variables = [lmi_variables x_lmi_variables(i)];
		j = j+1;
	end
end   
Y.dim(1) = size(trilX,1);
Y.dim(2) = size(trilX,2);
Y.lmi_variables = lmi_variables;
% Reset info about conic terms
Y.conicinfo = [0 0];
Y = clean(Y);