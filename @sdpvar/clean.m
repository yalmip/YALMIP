function Z=clean(X,tol)
%CLEAN Remove terms with small coefficients
%
% Z = clean(X,tol) removes all variables with a coefficient smaller than tol

X.extra.opname='';
if nargin==1
    temp = any(X.basis,1);
    temp = temp(2:end);
    index = find(temp);
else
    X.basis(abs(X.basis)<tol)=0;
    index = find(any(abs(X.basis(:,2:end))>tol,1));
end
if ~isempty(index)
    Z = X;
    if length(index)~=length(Z.lmi_variables)
        Z.basis = Z.basis(:,[1 1+index]);
        Z.lmi_variables = Z.lmi_variables(index);
    end
    if ~isreal(Z.basis) & nargin==2
        basis_real = real(Z.basis);
        basis_imag = imag(Z.basis);
        basis_real(abs(basis_real)<tol) = 0;
        basis_imag(abs(basis_imag)<tol) = 0;
        Z.basis = basis_real + sqrt(-1)*basis_imag;
    end        
else
    Z = full(reshape(X.basis(:,1),X.dim(1),X.dim(2)));
end