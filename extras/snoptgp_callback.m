function [F,G] = snoptgp_callback(x,mmodel)

persistent model
if nargin > 1
    model = mmodel;
    return
end

[f,df] = fmincon_fungp(x,model);
[g,geq,dg,dgeq] = fmincon_congp(x,model);
F = [f;g;geq];
n = length(x);
m = length(g);
p = length(geq);
G = [reshape(df,1,n);dg';dgeq'];