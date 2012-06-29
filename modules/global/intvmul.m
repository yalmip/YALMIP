function C = intvmul(A,B)
a = A(1);
b = A(2);
c = B(1);
d = B(2);
C = [min([a*c, a*d, b*c, b*d]), max([a*c, a*d, b*c, b*d])];
