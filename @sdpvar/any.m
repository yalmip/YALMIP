function anys = any(x)
%ANY (overloaded)

% Author Johan Löfberg 
% $Id: any.m,v 1.2 2006-07-26 20:17:57 joloef Exp $  

x_base = x.basis;
anys = full(sum(abs(x.basis),2)>0);
anys = reshape(anys,x.dim(1),x.dim(2));