function X = constraintifyUncertaintyDeclared(X)

try
    v = struct(X).extra.distribution.type;
    if isequal(v,'deterministic')
        X.typeflag = 15;
    else
        X.typeflag = 16;
    end
catch
end