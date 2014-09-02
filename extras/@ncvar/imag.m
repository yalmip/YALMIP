function X = imag(X)
%IMAG (overloaded)

X.basis = imag(X.basis);
X = clean(X);
if isa(X,'sdpvar')
   X.conicinfo = [0 0];
end