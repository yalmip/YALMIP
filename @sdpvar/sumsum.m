function Y=sumsum(Y)
%SUM (overloaded)

Y.basis = sum(Y.basis,1);
Y.dim(1) = 1;
Y.dim(2) = 1;
% Reset info about conic terms
Y.conicinfo = [0 0];
Y = clean(Y);
