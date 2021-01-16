function C = ttimes(A,B)
% TTIMES Tropical multiplication (max-plus algebra)
%
% C = TTIMES(A,B) Computes the tropical multiplication of the matrices A
% and B, AB(i,j) = max(A(i,:) + B(:,j)');
%
% See also tplus

n = size(A,1);
m = size(B,2);
C = reshape(max(kron(ones(m,1),A)+kron(B',ones(n,1)),[],2),n,m);
