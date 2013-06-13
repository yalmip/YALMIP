function y = times(X,Y)
%TIMES (overloaded)

% Check dimensions
[n,m]=size(X);
if ~((prod(size(X))==1) | (prod(size(Y))==1))
    if ~((n==size(Y,1) & (m ==size(Y,2))))
        error('Matrix dimensions must agree.')
    end
end

% Reshape if one argument is scalar
dX = size(X);
dY = size(Y);
if max(dX)==1 & max(dY)>1
    X = X*ones(dY);
    dX = size(X);
end
if max(dY)==1 & max(dX)>1
    Y = Y*ones(dX);
    dY = size(Y);
end

% Just loop and call mtimes instead
y = [];
X = X(:);
Y = Y(:);
for i = 1:length(X)
    y = [y;extsubsref(X,i,1)*extsubsref(Y,i,1)];
end
y = reshape(y,dX);
return
