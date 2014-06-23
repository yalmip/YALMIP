function [g,geq,dg,dgeq,xevaled] = fmincon_con(x,model,xevaled)

global latest_xevaled
global latest_x_xevaled
% Early bail for linear problems
g = [];
geq = [];
dg = [];
dgeq = [];
if model.linearconstraints
    xevaled = [];
    return
end

if nargin<3
    if isequal(x,latest_x_xevaled)
        xevaled = latest_xevaled;
    else
        xevaled = zeros(1,length(model.c));
        xevaled(model.linearindicies) = x;
        xevaled = apply_recursive_evaluation(model,xevaled);
        latest_x_xevaled = x;
        latest_xevaled = xevaled;
    end
end

if model.nonlinearinequalities
    g = full(model.Anonlinineq*xevaled(:)-model.bnonlinineq);
end

if model.nonlinearequalities
    geq = full(model.Anonlineq*xevaled(:)-model.bnonlineq);
end

dgAll_test = [];

if nargout == 2 || ~model.derivative_available
    return
elseif ~isempty(dgAll_test) & isempty(model.evalMap)
    dgAll = dgAll_test;
elseif isempty(model.evalMap) & (model.nonlinearinequalities | model.nonlinearequalities)
    n = length(model.c);
    linearindicies = model.linearindicies;    
    xevaled = zeros(1,n);
    xevaled(linearindicies) = x;
    % FIXME: This should be vectorized
    
    news = model.fastdiff.news;
    allDerivemt = model.fastdiff.allDerivemt;
    c = model.fastdiff.c;
    
    if model.fastdiff.univariateDifferentiates
        zzz = c.*(x(model.fastdiff.univariateDiffMonom).^model.fastdiff.univariateDiffPower);
    else
        %  X = repmat(x(:)',length(c),1);
        O = ones(length(c),length(x));
        nz = find(allDerivemt);
        %  O(nz) = X(nz).^allDerivemt(nz);
        O(nz) = x(ceil(nz/length(c))).^allDerivemt(nz);        
        zzz = c.*prod(O,2);               
    end
                          
    newdxx = model.fastdiff.newdxx;
    newdxx(model.fastdiff.linear_in_newdxx) = zzz;                
    newdxx = newdxx';
    
    if ~isempty(model.Anonlineq)    
        dgAll = model.Anonlineq*newdxx;
    else
        dgAll = [];
    end
    if ~isempty(model.Anonlinineq)    
        aux = model.Anonlinineq*newdxx;
        dgAll = [dgAll;aux];
    end    
    
else
    allA = [model.Anonlineq;model.Anonlinineq];
    requested = model.fastdiff.requested;
    dx = apply_recursive_differentiation(model,xevaled,requested,model.Crecursivederivativeprecompute);
    dgAll = allA*dx;
end

if  model.nonlinearequalities
    dgeq = dgAll(1:size(model.Anonlineq,1),:)';
end
if model.nonlinearinequalities
    dg = dgAll(size(model.Anonlineq,1)+1:end,:)';
end
