function y = addgkyp(X,Y)

sx = size(X);
sy = size(Y);

if ~((prod(sx)==1) | (prod(sy)==1))
    if ~all(size(X)==size(Y))
        error('Dimension mismatch');
    end
end

if isa(Y,'double') | (isa(Y,'sdpvar') & ~is(Y,'gkyp'))   
    X.extra.M =  X.extra.M+Y;
    X.typeflag = 0;
    y = X+Y;
    if isa(y,'double')
        return
    else
        y.typeflag = 40;
    end
else
    y = X;
    m = length(Y.extra.K);
    y.extra.M = y.extra.M + Y.extra.M;
    for i = 1:m
        y.extra.K{end+1} = Y.extra.K{i};
        y.extra.Phi{end+1} = Y.extra.Phi{i};
        y.extra.P{end+1} = Y.extra.P{i};
        y.extra.negated(end+1) = Y.extra.negated(i);
    end
    y.typeflag = 0;
    Y.typeflag = 0;
    y = y+Y;
    if isa(y,'double')
        return
    else
        y.typeflag = 40;
    end
end




