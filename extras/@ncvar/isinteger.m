function isint=isinteger(X)
%ISINTEGER Check if (part of) a variable is integer

% Author Johan Löfberg 
% $Id: isinteger.m,v 1.1 2006-08-10 18:00:20 joloef Exp $  

isint = any(ismember(getvariables(X),yalmip('intvariables')));
