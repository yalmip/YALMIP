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
    allA = [model.Anonlineq;model.Anonlinineq];
    dgAll = [];
    n = length(model.c);
    linearindicies = model.linearindicies;
    mtNonlinear = model.monomtable(model.nonlinearindicies,:);
    xevaled = zeros(1,n);
    xevaled(linearindicies) = x;
    X = repmat(xevaled,size(mtNonlinear,1),1);
    % FIXME: This should be vectorized
        
     news = model.fastdiff.news;
     allDerivemt = model.fastdiff.allDerivemt;
     c = model.fastdiff.c;
        
    zzz = c.*prod(repmat(x(:)',length(c),1).^allDerivemt,2);
    newdxx = spalloc(length(linearindicies),max(linearindicies),length(linearindicies));
    for i = 1:length(c)    
        newdxx(news(i,2),model.nonlinearindicies(news(i,1)))=zzz(i);%c(i)*prod(x(:)'.^allDerivemt(i,:));
    end
    for i = 1:length(linearindicies)
        newdxx(i,linearindicies(i)) = 1;
    end
    dgAll = allA*newdxx';
        
else
    allA = [model.Anonlineq;model.Anonlinineq];
    requested = any(allA',2);
    [i,j,k] = find((model.deppattern(find(requested),:)));
    requested(j) = 1;
    dx = apply_recursive_differentiation(model,xevaled,requested,model.Crecursivederivativeprecompute);    
    dgAll = allA*dx;
end

if  model.nonlinearequalities
    dgeq = dgAll(1:size(model.Anonlineq,1),:)';
end
if model.nonlinearinequalities
    dg = dgAll(size(model.Anonlineq,1)+1:end,:)';
end
