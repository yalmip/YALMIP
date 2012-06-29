function Y=sum(Y)
%SUM (overloaded)

% Author Johan Löfberg 
% $Id: sumsum.m,v 1.1 2006-08-10 18:00:22 joloef Exp $   

Y.basis = sum(Y.basis,1);
Y.dim(1) = 1;
Y.dim(2) = 1;
% Reset info about conic terms
Y.conicinfo = [0 0];
Y = clean(Y);
