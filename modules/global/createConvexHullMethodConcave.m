function [Ax,Ay,b,K] = createConvexHullMethodConcave(xL,xU,f,df,dummy)
% dummy is operator struct sent by blackbox
xM = (xL+xU)/2;
fL = real(f(xL));
fM = real(f(xM));
fU = real(f(xU));
dfL = real(df(xL));
dfM = real(df(xM));
dfU = real(df(xU));
x1 = xL;x2 = xU;
goal_derivative = (fU-fL)/(xU-xL);
for i = 1:4
    if dfM<goal_derivative
        x2 = xM;xM = (x1+x2)/2;fM = real(f(xM));dfM = real(df(xM));
    else
        x1 = xM;xM = (x1+x2)/2;fM = real(f(xM));dfM = real(df(xM));
    end
end
[Ax,Ay,b,K] = convexhullConcave(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);