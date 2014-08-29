function F = eq(X,Y)
% Internal class for constraint lists

superiorto('sdpvar');
superiorto('double');

if isa(X,'sdpvar')
    if is(X,'binary')
        F=iff(X,Y);
        return;
    end
end
if isa(Y,'sdpvar')
    if is(Y,'binary')
        F=iff(Y,X);
        return;
    end
end

error('Equalities can not be used in double-sided constraints')
