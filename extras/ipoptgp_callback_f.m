function f = ipoptgp_callback_f(x,prob)

global latest_x_f
global latest_df

x = x(:);
[f,latest_df] = fmincon_fungp(x,prob);
latest_x_f = x;