function Z = plus(X,Y)

if isa(X,'probability') && isa(Y,'probability')
    Z = X;
    for i = 1:length(Y.Weight)
        Z.Weight{end+1} = Y.Weight{i};
        Z.Risk{end+1} = Y.Risk{i};
        Z.Offset{end+1} = Y.Offset{i};
        Z.Constraint{end+1} = Y.Constraint{i};
    end
elseif isa(X,'probability') && (isa(Y,'sdpvar') || isa(Y,'double'))
    Z = X;
    Z.Weight{end+1} = 1;
    Z.Risk{end+1} = 0;
    Z.Offset{end+1} = Y;
    Z.Constraint{end+1} = [];
elseif isa(Y,'probability') && (isa(X,'sdpvar') || isa(X,'double'))
    Z = Y;
    Z.Weight{end+1} = 1;
    Z.Risk{end+1} = 0;
    Z.Offset{end+1} = X;
    Z.Constraint{end+1} = [];
end
