function Z = plus(X,Y)
% PLUS (overloaded)

dim = [];
if isa(X,'ndsdpvar')
    if isa(Y,'double') && numel(Y)==1
        % Quick return for scalar  + ndsdpvar
        Z = X;Z.basis(:,1) = Z.basis(:,1) + Y;
        Z.conicinfo = [0 0];
        Z.extra.opname='';
        return
    else
        dim = X.dim;
        X = sdpvar(X);
    end
elseif isa(X,'double')
    X = X(:);
end

if isa(Y,'ndsdpvar')
    if isa(X,'double') && numel(X)==1
        % Quick return for scalar  + ndsdpvar
        Z = Y;Z.basis(:,1) = Z.basis(:,1) + X;
        Z.conicinfo = [0 0];
        Z.extra.opname='';
        return
    else
        if isempty(dim)
            dim = Y.dim;
        else
            if isequal(dim,Y.dim)
            else
                error('Dimension mismatch in nD addition')
            end
        end
        Y = sdpvar(Y);
    end
elseif isa(Y,'double')
    Y = Y(:);
end

Z = X + Y;
Z = reshape(Z,dim);