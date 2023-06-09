function M = permutation(P)
% PERMUTATION
%
% M = permutation(X)
%
% Creates the model [sum(P,1)==1,sum(P,2)==1,binary(P)]

M = [sum(P,1)==1,sum(P,2)==1,binary(P)];