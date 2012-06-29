function Z = conj(Y)
%CONJ (overloaded)

% Author Johan Löfberg 
% $Id: conj.m,v 1.3 2006-01-26 13:44:13 joloef Exp $   

Z = Y;
Z.basis = conj(Z.basis);
% Reset info about conic terms
Z.conicinfo = [0 0];