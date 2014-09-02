function X = double(X)
% DOUBLE (overloaded)

X = reshape(double(sdpvar(X)),X.dim);