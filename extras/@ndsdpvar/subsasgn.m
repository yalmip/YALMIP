function X = subsasgn(X,I,Y)
% SUBSASGN (overloaded)

base = reshape(1:size(X.basis,1),X.dim);
base = subsref(base,I);

if isa(Y,'ndsdpvar')
    Y = sdpvar(Y);
elseif isa(Y,'double')
    Y = Y(:);
end

dim = X.dim;
X = sdpvar(X);
X(base) = Y;
X = reshape(X,dim);
X = clean(X);