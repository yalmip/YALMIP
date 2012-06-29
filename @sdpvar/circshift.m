function X = circshift(varargin)
% CIRCSHIFT (overloaded)

% Author Johan Löfberg
% $Id: circshift.m,v 1.2 2006-07-26 20:17:57 joloef Exp $

Y = varargin{1};
X = Y;
X.basis = [];
for i = 1:size(Y.basis,2)
    base = reshape(full(Y.basis(:,i)),[X.dim(1) X.dim(2)]);
    base = circshift(base,varargin{2:end});
    X.basis = [X.basis sparse(base(:))];
end
X.dim = size(base);
X.conicinfo = [0 0];
X = clean(X);