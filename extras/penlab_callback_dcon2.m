function [G,model] = penlab_callback_dcon2(x,k,model)

global latest_x_g
global latest_G
global latest_g
x = x(:);

G = spalloc(length(model.c),length(model.c),0);
