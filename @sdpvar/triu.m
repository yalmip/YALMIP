function X=triu(X,r)
%TRIU (overloaded)

% Author Johan Löfberg 
% $Id: triu.m,v 1.6 2006-07-26 20:17:58 joloef Exp $   


if nargin==1
    r = 0;
end

n = X.dim(1);
m = X.dim(2);

ind = reshape(1:n*m,n,m);
ind = triu(ind,r);
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

if nargin==1
    r = 0;
end

Y = X;
x_lmi_variables = X.lmi_variables;
lmi_variables = [];
n = X.dim(1);
m = X.dim(2);
triuX = triu(reshape(X.basis(:,1),n,m),r);
Y.basis = triuX(:);

j = 1;
for i = 1:length(x_lmi_variables)
	triuX = triu(reshape(X.basis(:,i+1),n,m),r);
	if (norm(triuX,inf)>0)
		Y.basis(:,j+1) = triuX(:);
		lmi_variables = [lmi_variables x_lmi_variables(i)];
		j = j+1;
	end
end   
Y.dim(1) = size(triuX,1);
Y.dim(2) = size(triuX,2);
Y.lmi_variables = lmi_variables;
% Reset info about conic terms
Y.conicinfo = [0 0];
Y = clean(Y);