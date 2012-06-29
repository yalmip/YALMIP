function X = sum(varargin)
% SUM (overloaded)

% Author Johan Löfberg
% $Id: sum.m,v 1.5 2006-11-07 08:18:28 joloef Exp $

Y = varargin{1};
X = Y;
X.basis = [];
for i = 1:size(Y.basis,2)
    base = reshape(full(Y.basis(:,i)),X.dim);
    base = sum(base,varargin{2:end});
    X.basis = [X.basis sparse(base(:))];
end
X.dim = size(base);
X.conicinfo = [0 0];
X = clean(X);