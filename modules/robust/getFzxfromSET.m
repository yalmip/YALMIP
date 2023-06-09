function [F0,Fz,Fx,Fxz] = getFzxfromSET(F,z,x)

F = sdpvar(F);
F0 = -getbasematrix(F,0);
zvars = getvariables(z);
xvars = getvariables(x);

for i = 1:length(xvars)
    Fx{i} = -getbasematrix(F,xvars(i));
end

for i = 1:length(zvars)
    Fz{i} = -getbasematrix(F,zvars(i));
end

for i = 1:length(xvars)   
    for j = 1:length(zvars)
        Fxz{i,j} = -getbasematrix(F,getvariables(x(i)*z(j)));
    end
end
