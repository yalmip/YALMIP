function c = unionstripped(a,b)
%UNIONSTRIPPED  Internal function (version without checkings etc.)

% Author Johan Löfberg
% $Id: unionstripped.m,v 1.1 2004-11-24 09:13:05 johanl Exp $

c = uniquestripped([a(:)' b(:)']);


