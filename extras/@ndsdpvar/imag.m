function X = imag(X)
% REAL (overloaded)

X.basis = imag(X.basis);
X.conicinfo = [0 0];
X = clean(X);