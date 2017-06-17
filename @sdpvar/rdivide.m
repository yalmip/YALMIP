function y = rdivide(X,Y)
%RDIVIDE (overloaded)

% Check dimensions
[nx,mx] = size(X);
[ny,my] = size(Y);
if ~((prod(size(X))==1) | (prod(size(Y))==1))
    if ~((nx==ny & (mx ==my)))
        error('Matrix dimensions must agree.')
    end
end

% Quick exit for simple case X/scalar
if isnumeric(Y) & prod(size(Y))==1
    y = X;
    y.basis = y.basis/Y;
    % Reset info about conic terms
    y.conicinfo = [0 0];
    return
end

if isa(X,'sdpvar') & isnumeric(Y)
    y = X.*(1./Y);
    return
end

% normalize scalar./matrix and matrix./scalar
[nx,mx] = size(X);
[ny,my] = size(Y);
if prod(size(X)) == 1 & prod(size(Y))~=1
    X = repmat(X,ny,my);
end;
if prod(size(Y)) == 1 & prod(size(X))~=1
    Y = repmat(Y,nx,mx);
end

% power is optimized already for simple cases
y = X.*(Y.^(-1));

% Reset info about conic terms
if isa(y,'sdpvar')
    y.conicinfo = [0 0];
end