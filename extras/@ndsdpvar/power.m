function X = power(X,d)
% POWER (overloaded)

% Vectorize x if d is vector
if numel(X)==1 & (numel(d)>1)
    X = X.*ones(size(d));
end
% Vectorize if x is a vector
if numel(d)==1 & (numel(X)>1)
    d = d.*ones(size(X));
end
s = size(X);
if isa(X,'sdpvar')
    X = sdpvar(X);
else
    X = X(:);
end
d = reshape(d,[],1);
X = power(X,d);
X = reshape(X,s);