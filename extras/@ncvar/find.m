function indicies = find(x)
%FIND (overloaded)

% Author Johan Löfberg 
% $Id: find.m,v 1.1 2006-08-10 18:00:20 joloef Exp $   

base = x.basis;
vars = x.lmi_variables;
indicies = find(any(base,2));
