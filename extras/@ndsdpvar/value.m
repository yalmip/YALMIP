function X = value(X)
% VALUE (overloaded)

X = reshape(value(sdpvar(X)),X.dim);