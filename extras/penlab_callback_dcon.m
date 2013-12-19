function [G,model] = penlab_callback_dcon(x,model)

global latest_x_g
global latest_G
global latest_g
x = x(:);

G = [model.Aeq;model.A];
