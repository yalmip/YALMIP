function Z = uplus(Y)
%UPLUS (overloaded)

% Author Johan Löfberg 
% $Id: uplus.m,v 1.1 2006-08-10 18:00:23 joloef Exp $   

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

