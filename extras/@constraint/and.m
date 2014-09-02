function F = and(X,Y)
% Internal class for constraint list

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

F = X & Y;