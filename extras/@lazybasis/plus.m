function y = plus(X,Y)
%PLUS (overloaded)

% Author Johan Löfberg
% $Id: plus.m,v 1.1 2005-10-12 16:05:54 joloef Exp $

if isa(X,'lazybasis');
    X = double(X);
end
if isa(Y,'lazybasis');
    Y = double(Y);
end
y = X + Y;

