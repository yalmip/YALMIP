function  out = isequal(X,Y,check)
%ISEQUAL (overloaded)

% Author Johan Löfberg 
% $Id: isequal.m,v 1.3 2007-04-26 12:52:22 joloef Exp $   

if nargin == 3
    if (isa(X,'sdpvar') & isa(Y,'sdpvar'))
        out = isequal(X.basis,Y.basis) & isequal(X.lmi_variables,Y.lmi_variables);
    else
        out = 0;
    end
else
    if (isa(X,'sdpvar') & isa(Y,'sdpvar'))
        out = isequal(struct(X),struct(Y));
    else
        out = 0;
    end
end
