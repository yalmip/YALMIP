function x = double(x);
%DOUBLE (Overloaded)

% Author Johan Löfberg
% $Id: double.m,v 1.1 2005-04-29 16:28:03 joloef Exp $

x = double(x.Evaluated{1});

if issymmetric(x)
    x = min(eig(x)) >= 0;
else
    x = min(x(:)) >= 0;
end
