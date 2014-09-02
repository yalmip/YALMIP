function  out = isequal(X,Y)
%ISEQUAL (overloaded)

if (isa(X,'ncvar') & isa(Y,'ncvar'))
    out = isequal(struct(X),struct(Y));
else
	out = false;
end
	