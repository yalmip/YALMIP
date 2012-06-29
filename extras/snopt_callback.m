function [F,G] = snopt_callback(x,mmodel)

persistent model
if nargin > 1
    model = mmodel;
    return
end

if model.SimpleLinearObjective | model.SimpleQuadraticObjective
    [f,df] = fmincon_fun(x,model);
    [g,geq,dg,dgeq] = fmincon_con(x,model);
else
    [f,df,xevaled] = fmincon_fun(x,model);
    [g,geq,dg,dgeq] = fmincon_con(x,model,xevaled);
end

F = [f;g;geq];
n = length(x);
m = length(g);
p = length(geq);
G = [reshape(df,1,n);dg';dgeq'];
if ~isempty(model.A)
    F = [F;model.A*x - model.b];
    G = [G;model.A];
end
if ~isempty(model.Aeq)
    F = [F;model.Aeq*x - model.beq];
    G = [G;model.Aeq];
end