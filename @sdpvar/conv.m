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
ii = [];
jj = [];
ss = [];
for i = 1:length(x_lmi_variables)+1
    x=reshape(X.basis(:,i),n,m);
    z=conv(full(x),full(Y));
    [iiz,jjz,ssz] = find(z(:));
    ii = [ii iiz(:)'];
    jj = [jj jjz(:)'+i-1];
    ss = [ss ssz(:)'];    
end  
Z.basis = sparse(ii,jj,ss);
Z.dim(1) = size(z,1);
Z.dim(2) = size(z,2);
Z = clean(Z);
% Reset info about conic terms
if isa(Z,'sdpvar')
    Z.conicinfo = [0 0];
    Z = flush(Z);
end