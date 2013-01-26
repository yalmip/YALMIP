function X=brutepersp(X,t,y)

X.basis = [zeros(size(X.basis,1),1) X.basis];
X.lmi_variables = [t y];
[i,j] = sort(X.lmi_variables);
X.basis = X.basis(:,[1 1+j]);
X.lmi_variables = i;
X.conicinfo = [0 0];