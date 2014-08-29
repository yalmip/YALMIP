function X = dual2cell(dual_vec,K)
%DUAL2CELL Internal function for organizing dual data 

row = 1;

X.f = dual_vec(row:row+K.f-1);
row = row + K.f;

X.l = dual_vec(row:row+K.l-1);
row = row + K.l;

for k = 1:length(K.q)
    X.q{k} = dual_vec(row:row+K.q(k)-1);
    row = row + K.q(k);
end

for k = 1:length(K.r)
    X.r{k} = dual_vec(row:row+K.r(k)-1);
    row = row + K.r(k);
end

for k = 1:length(K.s)
    X.s{k} = mat(dual_vec(row:row+K.s(k)^2-1));
    row = row + K.s(k)^2;
end

function Y = mat(X)
n = sqrt(length(X));
Y = reshape(X,n,n);

