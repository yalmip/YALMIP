function Z=conv(X,Y)
%CONV (overloaded)

if isnumeric(X)
    temp = X;
    X = Y;
    Y = temp;
end
if isa(X,'sdpvar') & isa(Y,'sdpvar')
    error('nonlinear CONV not supported yet. Make feature request');
end
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
Z = clean(Z);
% Reset info about conic terms
if isa(Z,'sdpvar')
    Z.conicinfo = [0 0];
    Z = flush(Z);
end