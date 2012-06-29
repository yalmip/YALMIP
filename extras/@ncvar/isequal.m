function  out = isequal(X,Y)
%ISEQUAL (overloaded)

% Author Johan Löfberg 
% $Id: isequal.m,v 1.1 2006-08-10 18:00:20 joloef Exp $   

if (isa(X,'ncvar') & isa(Y,'ncvar'))
    out = isequal(struct(X),struct(Y));
else
	out = 0;
end
	