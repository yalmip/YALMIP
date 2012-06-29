function y = getbasematrix(x,ind)
%GETBASEMATRIX (overloaded sdpvar/getbasematrix on double)

% Author Johan Löfberg
% $Id: getbasematrix.m,v 1.2 2004-07-02 08:17:30 johanl Exp $

if ind == 0
    y = x;
else
    y = spalloc(size(x,1),size(x,2),0);
end
