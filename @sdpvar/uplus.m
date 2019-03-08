function Z = uplus(Y)
%UPLUS (overloaded)

disp('Most likely you meant to write a + b, but you wrote a +b')
disp('This can easily lead to bugs, as [a +b] is a vector with two elements')
disp('If you really want to use unitary plus, you will have to edit sdpvar/uplus')
disp('and delete this message')
if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

switch Y.typeflag
case {0,5}
	Z = Y;
case {1,3,4}
	error('Relational objects canot be manipulated')
otherwise
	error('Please report internal bug in sdpvar/uplus')    
end

