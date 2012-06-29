function Y=sum(Y)
%SUM (overloaded)

% Author Johan Löfberg 
% $Id: sumsum.m,v 1.3 2006-07-26 20:17:58 joloef Exp $   

Y.basis = sum(Y.basis,1);
Y.dim(1) = 1;
Y.dim(2) = 1;
% Reset info about conic terms
Y.conicinfo = [0 0];
Y = clean(Y);
