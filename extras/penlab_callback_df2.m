function [H,model] = penlab_callback_df2(x,model)

global latest_x_f
global latest_df

H = speye(2);