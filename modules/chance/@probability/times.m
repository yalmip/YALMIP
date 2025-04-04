function Z = times(X,Y)

if isa(X,'probability') && isa(Y,'probability')
    error()
elseif isa(X,'probability') && (isa(Y,'sdpvar') || isa(Y,'double'))
    Z = X;
    if isa(Y,'double') && (Y<0)
        error('Can only multiply probability objects with positive constants');
    end
    for i = 1:length(X.Weight)
        Z.Weight{i} = Z.Weight{i}*Y;
    end
elseif isa(Y,'probability') && (isa(X,'sdpvar') || isa(X,'double'))
    Z = Y;
    if isa(X,'double') && (X<0)
        error('Can only multiply probability objects with positive constants');
    end
    for i = 1:length(Y.Weight)
        Z.Weight{i} = Z.Weight{i}*X;
    end
end
