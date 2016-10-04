function [g,geq,dg,dgeq,xevaled] = fmincon_con(x,model,xevaled)

if isempty(model.lift)    
    [g,geq,dg,dgeq,xevaled] = fmincon_con(x,model);
else
    xlift = zeros(length(model.linearindicies),1);
    xlift(model.lift.linearIndex) = x;
    xlift(model.lift.liftedIndex) = model.lift.d + model.lift.T*x;
    
    % Call the computational kernel which works in the fully expanded
    % normalized format
    [g,geq,dg,dgeq,xevaled] = fmincon_con(xlift,model);
     
  
    [f,df,xevaledout] = fmincon_fun(xlift,model);
    %Now map gradient to exposed varaibles to fmincon
    if ~isempty(dg)
        dg = dg(model.lift.linearIndex,:) + model.lift.T'*dg(model.lift.liftedIndex,:);
    end
    if ~isempty(dgeq)
        dgeq = dgeq(model.lift.linearIndex,:) + model.lift.T'*dgeq(model.lift.liftedIndex,:);
    end
end
