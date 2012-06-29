function f = ipopt_callback_f(x,model_in)

global latest_x
global latest_df

persistent model
if nargin>1
    model = model_in;
    return
end

% We compute both function value and gradient. If pennon calls later to get
% gradient, we have saved in a global variable (yea, I know, globals suck)
x = x(:);
[f,latest_df] = fmincon_fun(x,model);
latest_x = x;
