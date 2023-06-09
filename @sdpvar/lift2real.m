function X = lift2real(X)

reF=real(X);
imF=imag(X);
X = [reF -imF;imF reF];