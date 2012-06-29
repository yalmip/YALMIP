function y = mpower(x,d)
%MPOWER (overloaded)

% Author Johan Löfberg
% $Id: mpower.m,v 1.2 2006-08-11 11:48:15 joloef Exp $

%Sanity check
if prod(size(d))>1 | ~((fix(d) == d)) | ~(d>=0)
    error('The power must be scalar.');
end
if x.dim(1)~=x.dim(2) 
    error('Matrix must be square.')
end
    
switch d
    case 0
            y = eye(x.dim(1),x.dim(2))^0;
    case 1
        y = x;
    case 2
        y = x * x;
    otherwise
        y = x*mpower(x,d-1);
end
