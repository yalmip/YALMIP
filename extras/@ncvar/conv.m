function Z=conv(X,Y)
%CONV (overloaded)

% Author Johan Löfberg 
% $Id: conv.m,v 1.1 2006-08-10 18:00:19 joloef Exp $   


x_lmi_variables = X.lmi_variables;
n = X.dim(1);
m = X.dim(2);
Z=X;
Z.basis=[];
for i = 1:length(x_lmi_variables)+1
    x=reshape(X.basis(:,i),n,m);
    z=conv(full(x),full(Y));
    Z.basis=[Z.basis z(:)];
end   
Z.dim(1) = size(z,1);
Z.dim(2) = size(z,2);
% Reset info about conic terms
Z.conicinfo = [0 0];