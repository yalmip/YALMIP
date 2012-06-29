function Z = plus(X,Y)
% PLUS (overloaded)

% Author Johan Löfberg
% $Id: plus.m,v 1.2 2006-07-13 19:40:59 joloef Exp $

dim = [];
if isa(X,'ndsdpvar')
    dim = X.dim;
    X = sdpvar(X);
elseif isa(X,'double')
    X = X(:);
end

if isa(Y,'ndsdpvar')
    if isempty(dim)
        dim = Y.dim;
    else
        if isequal(dim,Y.dim)
        else
            error('Dimension mismatch in nD addition')
        end
    end
    Y = sdpvar(Y);
elseif isa(Y,'double')
    Y = Y(:);
end

Z = X + Y;
Z = reshape(Z,dim);