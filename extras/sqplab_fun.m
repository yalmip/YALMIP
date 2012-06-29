function [f,df] = sqplab_fun(x)

global SQPLABDATA

xevaled = zeros(1,length(SQPLABDATA.interfacedata.c));
xevaled(SQPLABDATA.linearindicies) = x;

% Experimental support for arbitrary functions
if SQPLABDATA.evalinobjective
    if SQPLABDATA.monominobjective
        if ~isempty(SQPLABDATA.bilinears)
            xevaled(SQPLABDATA.bilinears(:,1)) = xevaled(SQPLABDATA.bilinears(:,2)).*xevaled(SQPLABDATA.bilinears(:,3));
        else
            xevaled(SQPLABDATA.nonlinearindicies) = prod(repmat(xevaled,length(SQPLABDATA.nonlinearindicies),1).^SQPLABDATA.monomtable(SQPLABDATA.nonlinearindicies,:),2);
        end
    end
    for i = 1:length(SQPLABDATA.interfacedata.evalMap)
        arguments = SQPLABDATA.evalMap{i}.prearg;
        arguments{2} = xevaled(SQPLABDATA.interfacedata.evalMap{i}.variableIndex);
        xevaled(SQPLABDATA.interfacedata.evalVariables(i)) = feval(arguments{:});
    end
end

if SQPLABDATA.monominobjective
    if ~isempty(SQPLABDATA.bilinears)
        xevaled(SQPLABDATA.bilinears(:,1)) = xevaled(SQPLABDATA.bilinears(:,2)).*xevaled(SQPLABDATA.bilinears(:,3));
    else
        xevaled(SQPLABDATA.nonlinearindicies) = prod(repmat(xevaled,length(SQPLABDATA.nonlinearindicies),1).^SQPLABDATA.monomtable(SQPLABDATA.nonlinearindicies,:),2);
    end
end

xevaled = xevaled(:);
if SQPLABDATA.SimpleLinearObjective
    f = SQPLABDATA.interfacedata.c'*xevaled;
else
    f = (SQPLABDATA.interfacedata.c'+xevaled'*SQPLABDATA.interfacedata.Q)*xevaled;
end

if SQPLABDATA.SimpleLinearObjective
    df = SQPLABDATA.interfacedata.c(SQPLABDATA.linearindicies);
elseif SQPLABDATA.SimpleQuadraticObjective
    df = SQPLABDATA.interfacedata.c(SQPLABDATA.linearindicies) + 2*SQPLABDATA.interfacedata.Q(SQPLABDATA.linearindicies,SQPLABDATA.linearindicies)*x;
elseif SQPLABDATA.SimpleNonlinearObjective
    df = [];
    for i = 1:length(SQPLABDATA.linearindicies)
        xevaled = zeros(1,length(SQPLABDATA.interfacedata.c));
        xevaled(SQPLABDATA.linearindicies) = x;
        mt = SQPLABDATA.monomtable;
        oldpower = mt(:,SQPLABDATA.linearindicies(i));
        mt(:,SQPLABDATA.linearindicies(i)) = mt(:,SQPLABDATA.linearindicies(i))-1;
        xevaled = prod(repmat(xevaled,size(mt,1),1).^mt,2);
        xevaled = xevaled(:)'.*oldpower';xevaled(isnan(xevaled))=0;
        df = [df;SQPLABDATA.interfacedata.c'*xevaled'];
    end
    df = df + 2*SQPLABDATA.interfacedata.Q(SQPLABDATA.linearindicies,SQPLABDATA.linearindicies)*x;
else
    df = [];
end
