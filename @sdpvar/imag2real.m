function Y = imag2real(Y)
%KRON (overloaded)

% Author Johan Löfberg
% $Id: imag2real.m,v 1.1 2007-03-02 09:51:35 joloef Exp $

if isreal(Y)
    y = Y;
    return
end

realBase = real(Y.basis);
imagBase = imag(Y.basis);

lmi_variables = getvariables(Y);
nv = length(lmi_variables);

% [re im;-im re] = kron(I,re) + kron([0 1;-1 0],im)
sparse_X1 = [1 0;0 1];
sparse_X2 = [0 1;-1 0];
temp = kron(sparse_X1,reshape(realBase(:,1),Y.dim(1),Y.dim(1)))+ kron(sparse_X2,reshape(imagBase(:,1),Y.dim(1),Y.dim(1)));
temp = temp(:);

Y.basis = temp(:);
for i = 1:nv
    temp1 = kron(sparse_X1,reshape(realBase(:,i+1),Y.dim(1),Y.dim(1)));
    temp2 = kron(sparse_X2,reshape(imagBase(:,i+1),Y.dim(1),Y.dim(1)));
    Y.basis(:,i+1) = temp1(:) + temp2(:);
end;

Y.dim(1) = size(temp1,1);
Y.dim(2) = size(temp1,2);
Y = clean(Y);
% Reset info about conic terms
if isa(Y,'sdpvar')
    Y.conicinfo = [0 0];
end