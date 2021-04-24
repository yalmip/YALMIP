function x = RecoverNonlinearSolverSolution(model,xout);

if isempty(model.nonlinearindicies)
    x = xout(:);
else
    x = zeros(length(model.c),1);
    for i = 1:length(model.linearindicies)
        x(model.linearindicies(i)) = xout(i);
    end
    x = x(1:length(model.c));
end
