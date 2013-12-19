function [df,model] = penlab_callback_df(x,model)

global latest_x_f
global latest_df

x = x(:);
if isequal(x,latest_x_f)
    df = latest_df;
else
    [f,df] = fmincon_fun(x,model);
end
df = df(:);