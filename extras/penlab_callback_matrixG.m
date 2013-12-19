function [G,model] = penlab_callback_matrixG(x,k,model)

global latest_xevaled
global latest_x_xevaled
x = x(:);

if isequal(x,latest_x_xevaled)
    xevaled = latest_xevaled;
else
    xevaled = zeros(1,length(model.c));
    xevaled(model.linearindicies) = x;
    xevaled = apply_recursive_evaluation(model,xevaled);
    latest_x_xevaled = x;
    latest_xevaled = xevaled;
end
    
vecG = model.vecF{k}*[1;xevaled];
G = reshape(vecG,model.K.s(k),model.K.s(k));
