function X = squeeze(X)
% squeeze (overloaded)

% Author Johan Löfberg
% $Id: squeeze.m,v 1.2 2006-07-26 09:09:07 joloef Exp $


dummy = reshape(ones(prod(X.dim),1),X.dim);
dummy = squeeze(dummy);
X.dim = size(dummy);
X = clean(X);
