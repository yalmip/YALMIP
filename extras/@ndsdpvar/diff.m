function X = diff(varargin)
% DIFF (overloaded)

% Author Johan Löfberg
% $Id: diff.m,v 1.2 2006-07-13 19:40:59 joloef Exp $

Y = varargin{1};
X = Y;
X.basis = [];
for i = 1:size(Y.basis,2)
    base = reshape(full(Y.basis(:,i)),X.dim);
    base = diff(base,varargin{2:end});
    X.basis = [X.basis sparse(base(:))];
end
X.dim = size(base);
