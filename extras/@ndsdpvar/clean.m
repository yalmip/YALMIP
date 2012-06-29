function Y=clean(Y)
%clean (overloaded)

% Author Johan Löfberg
% $Id: clean.m,v 1.1 2006-07-13 19:40:59 joloef Exp $

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
    Y = reshape(sdpvar(Y),Y.dim(1),Y.dim(2));
    Y = clean(Y);
else
    temp = any(Y.basis,1);
    temp = temp(2:end);
    index = find(temp);
    if ~isempty(index)
        Z = Y;
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
        Y = full(reshape(Y.basis(:,1),Y.dim));
    end
end