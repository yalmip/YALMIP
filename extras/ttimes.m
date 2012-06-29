function C = mptimes(A,B)
% TTIMES Tropical multiplication (max-plus algebra)
%
% C = TTIMES(A,B) Computes the tropical multiplication of the matrices A
% and B, AB(i,j) = max(A(i,:) + B(:,j)');
%
% See also tplus

% Author Johan Löfberg 
% $Id: ttimes.m,v 1.1 2006-11-24 13:39:43 joloef Exp $  

n = size(A,1);
m = size(B,2);
C = reshape(max(kron(ones(m,1),A)+kron(B',ones(n,1)),[],2),n,m);
% 
% C = [];
% for i = 1:size(A,1)
%     temp = [];
%     for j = 1:size(B,2)
%         temp = [temp max(A(i,:) + B(:,j)')];
%     end
%     C = [C;temp];
% end
