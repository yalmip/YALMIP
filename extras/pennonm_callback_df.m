function [nnz,ind,val] = pennonm_callback_df(x,model_in)

global latest_x
global latest_df

persistent model
if nargin>1
    model = model_in;
    return
end

x = x(:);
if isequal(x,latest_x)
    df = latest_df;
else
    [f,df] = fmincon_fun(x,model);
end

nnz = length(df);
ind = 1:length(df);
val = full(df);
