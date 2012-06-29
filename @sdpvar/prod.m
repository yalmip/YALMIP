function product=prod(X)
%PROD (overloaded)

% Author Johan L÷fberg
% $Id: prod.m,v 1.8 2008-01-16 22:14:12 joloef Exp $

% Ugly hack to fix a special user case prod([1e-7 1e-7 1e7 1e7].*[x1 x2 x3
% x4])
scale = 1;
if min(size(X))==1
    base = getbase(X);
    for i = 1:size(X.basis,1)
        si = norm(full(base(i,:)));
        if si<1
            scale = scale * si;
            base(i,:) = base(i,:)/si;
        end
    end
    X.basis = base;
end
product = 1;
for i = 1:length(X)
    pick = cell(1,1);pick{1}={i};
    product = product*subsref(X,struct('type','()','subs',pick));
end
product.basis = product.basis*scale;