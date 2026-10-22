function product=prod(X,dimen)
%PROD (overloaded)

if nargin == 2
    if dimen > 2 || size(X,dimen) == 1
        product = X;
        return
    end
    if dimen == 2
        product = prod(X')';
        return
    end
end

% Ugly hack to fix a special user case 
% such as prod([1e-7 1e-7 1e7 1e7].*[x1 x2 x3 x4])
scale = 1;
if min(size(X))==1
    base = getbase(X);
    for i = 1:size(X.basis,1)
        si = norm(full(base(i,:)));
        if si == 0
            product = 0;
            return
        elseif si<1
            scale = scale * si;
            base(i,:) = base(i,:)/si;
        end
    end
    X.basis = base;
end

if X.dim(1)==1 | X.dim(2)==1
      product = 1;    
    for i = 1:length(X)
        pick = cell(1,1);pick{1}={i};
        product = product*subsref(X,struct('type','()','subs',pick));
    end
    product = product*scale;
else
    vecProd = [];
    top = 1;
    for j = 1:X.dim(2)
        product = 1;
        for i = 1:X.dim(1)
            pick = cell(1,1);pick{1}={top};
            product = product*subsref(X,struct('type','()','subs',pick));
            top = top + 1;
        end
        vecProd = [vecProd;product];
    end
    product = reshape(vecProd,1,X.dim(2))*scale;
end
