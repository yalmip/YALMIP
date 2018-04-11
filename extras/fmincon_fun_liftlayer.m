function [f,df,xevaledout] = fmincon_fun_liftlayer(x,model)

if isempty(model.lift)    
    [f,df,xevaledout,dx] = fmincon_fun(x,model);
else
    xlift = zeros(length(model.linearindicies),1);
    xlift(model.lift.linearIndex) = x;
    xlift(model.lift.liftedIndex) = model.lift.d + model.lift.T*x;
    % Call the computational kernel which works in the fully expanded
    % normalized format.
    [f,df,xevaledout,dx] = fmincon_fun(xlift,model);
    %Now map gradient to exposed variables to fmincon
    df = df(:)';
    if ~isempty(df)
        df = df(model.lift.linearIndex) + df(model.lift.liftedIndex)*model.lift.T;
    end
end
