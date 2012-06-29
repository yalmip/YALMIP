function X=brutepersp(X,t,y)

% Author Johan Löfberg 
% $Id: brutepersp.m,v 1.3 2006-03-15 22:45:58 joloef Exp $   

X.basis = [zeros(size(X.basis,1),1) X.basis];
X.lmi_variables = [t y];
[i,j] = sort(X.lmi_variables);
X.basis = X.basis(:,[1 1+j]);
X.lmi_variables = i;
X.conicinfo = [0 0];