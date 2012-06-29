function sys=double(X)
%DOUBLE Returns current numerical value 

% Author Johan Löfberg 
% $Id: double.m,v 1.1 2005-10-12 16:05:54 joloef Exp $  

sys = sparse(X.iX,X.jX,X.sX,X.n,X.m);
%X.basis;