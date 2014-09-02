function isint=isinteger(X)
%ISINTEGER Check if (part of) a variable is integer

isint = any(ismember(getvariables(X),yalmip('intvariables')));
