function [f,err] = fmincon_funn(cbData,nRow,x,njdiff,dXjbase,reserved,inParam)

persistent params hash xevaled oldh F
if nargin == 7
    params = inParam;
    hash = randn(1,length(params.linearindicies));
    oldh = randn(1);
    xevaled = zeros(1,length(params.interfacedata.c));
    F = [];
    return
end

x = x(1:length(params.linearindicies));

if x(:)'*hash(:) ~= oldh
    oldh = x(:)'*hash(:);
    xevaled = zeros(1,length(params.interfacedata.c));
    xevaled(params.linearindicies) = x;

    % Experimental support for arbitrary functions
    if ~isempty(params.interfacedata.evalMap)
        
        pp = params.monomtable(params.nonlinearindicies,:);
        xevaled(params.nonlinearindicies) = prod(repmat(xevaled,length(params.nonlinearindicies),1).^pp,2);
        
        for i = 1:length(params.interfacedata.evalMap)
            arguments = {params.interfacedata.evalMap{i}.fcn,xevaled(params.interfacedata.evalMap{i}.variableIndex)};
            arguments = {arguments{:},params.interfacedata.evalMap{i}.arg{2:end-1}};
            xevaled(params.interfacedata.evalVariables(i)) = feval(arguments{:});
        end
    end
    xevaled(params.nonlinearindicies) = prod(repmat(xevaled,length(params.nonlinearindicies),1).^params.monomtable(params.nonlinearindicies,:),2);
    xevaled = xevaled(:);    
    F =  -params.F_struc*[1;xevaled];
end

if nRow == -1
    f = params.interfacedata.c'*xevaled+xevaled'*params.interfacedata.Q*xevaled;
else
    f = F(nRow + 1);%-params.F_struc(nRow + 1,:)*[1;xevaled];
end
err = 0;