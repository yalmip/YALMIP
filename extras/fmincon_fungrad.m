function df = fmincon_fun(x,params)

xevaled = zeros(1,length(params.interfacedata.c));
xevaled(params.linearindicies) = x;

% Experimental support for arbitrary functions
if ~isempty(params.interfacedata.evalMap)
    for i = 1:length(params.interfacedata.evalMap)
        xevaled(params.interfacedata.evalVariables(i)) = feval( params.interfacedata.evalMap{i}.fcn,xevaled(params.interfacedata.evalMap{i}.variableIndex));
    end
end

if nnz(params.interfacedata.c(params.nonlinearindicies)) == 0 & isempty(params.interfacedata.evalMap)
    %At most quadratic!
    df = params.interfacedata.c(params.linearindicies) + 2*params.interfacedata.Q(params.linearindicies,params.linearindicies)*x;
else
    error('not implmented')
end
    