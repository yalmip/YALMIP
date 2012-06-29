function isint=isinteger(X)
%ISINTEGER Check if (part of) a variable is integer

% Author Johan Löfberg 
% $Id: isinteger.m,v 1.2 2004-07-01 11:17:11 johanl Exp $  

isint = any(ismember(getvariables(X),yalmip('intvariables')));
