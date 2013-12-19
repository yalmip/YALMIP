function [G,model] = penlab_callback_con(x,model)

global latest_x_g
global latest_G
global latest_g
x = x(:);

G = [model.Aeq;model.A]*x;
