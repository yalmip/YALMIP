function [A,b] = randpolytope(n,m)

A = randn(n,m);
b = m*rand(n,1)-A*randn(m,1);