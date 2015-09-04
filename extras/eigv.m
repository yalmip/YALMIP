function [v,Model] = eigv(X)

n = size(X,1);
v = sdpvar(n,1);
Model = [trace(X)==sum(v), v(1:end-1) >= v(2:end)];
for i = 1:n-1 
    Model = [Model, sumk_generator(X,i,sum(v(1:i)))];
end
