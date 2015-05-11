function Y=clean(Y)
%clean (overloaded)

i = length(Y.dim);
while i>2
    if Y.dim(end) == 1
        Y.dim(end) = [];
        i = i - 1;
    else
        break
    end
end
if length(Y.dim) == 2;    
    Y = sdpvar([],[],[],[],[],[],[],struct(Y));
    Y = clean(Y);
else
    temp = any(Y.basis,1);
 %   temp = temp(2:end);
    index = find(temp);
    if ~isequal(index,1)%~isempty(index)
        Z = Y;
        if index(1)==1
            index = index(2:end)-1;
        else
            index = index-1;
        end
        if length(index)~=length(Y.lmi_variables)
            Y.basis = Y.basis(:,[1 1+index]);
            Y.lmi_variables = Y.lmi_variables(index);
        end
        if ~isreal(Y.basis) & nargin==2
            basis_real = real(Y.basis);
            basis_imag = imag(Y.basis);
            basis_real(abs(basis_real)<tol) = 0;
            basis_imag(abs(basis_imag)<tol) = 0;
            Y.basis = basis_real + sqrt(-1)*basis_imag;
        end
    else
        Y = reshape(full(Y.basis(:,1)),Y.dim);
    end
end