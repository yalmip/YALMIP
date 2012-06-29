function Z = conj(Y)
%CONJ (overloaded)

% Author Johan Löfberg 
% $Id: conj.m,v 1.1 2006-08-10 18:00:19 joloef Exp $   

Z = Y;
Z.basis = conj(Z.basis);
% Reset info about conic terms
Z.conicinfo = [0 0];