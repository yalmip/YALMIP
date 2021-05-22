function [f,df] = knitro_callback_f(x,model)

global latest_x_f
global latest_df
global latest_f

x = x(:);
if isequal(x,latest_x_f)
    df = latest_df;
    f = latest_f;
else
    [f,df] = fmincon_fun_liftlayer(x,model);
    latest_x_f = x;
    latest_f = f;
    latest_df = df;
end