function X = real(X)
% REAL (overloaded)

X.basis = real(X.basis);
X.conicinfo = [0 0];
X = clean(X);