function sys = and(X,Y)
%AND Overloaded
%   
%   See also   LMI

% TODO : Check if binaries etc
if isa(X,'sdpvar')
    X = true(X);
end
if isa(Y,'sdpvar')
    Y = true(Y);
end

if isa(X,'constraint')
    X = lmi(X);
end
if isa(Y,'constraint')
    Y = lmi(Y);
end

sys = plus(X,Y);