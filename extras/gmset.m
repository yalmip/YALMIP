function F = gmset(t,x1,x2)
%GMSET  Internal function used for MAXDET formulation

F =  (cone([t;(x1-x2)/2],(x1+x2)/2));