function b = uniquestripped(a)
%UNIQUESTRIPPED  Internal function (version without checkings etc.)

% Author Johan Löfberg
% $Id: uniquestripped.m,v 1.1 2004-11-24 09:13:05 johanl Exp $

b = sort(a(:)');
b = b(diff([b NaN])~=0);

