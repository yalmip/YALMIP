function [f,df,xevaledout,dx] = fmincon_fun(x,model)

global latest_xevaled
global latest_x_xevaled

% Apply the precomputed evaluation scheme (if necessary)
xevaled = zeros(1,length(model.c));
xevaled(model.linearindicies) = x;
if ~model.SimpleLinearObjective
    if isequal(x,latest_x_xevaled)
        xevaled = latest_xevaled;
    else       
        xevaled = apply_recursive_evaluation(model,xevaled);
        latest_xevaled = xevaled;
        latest_x_xevaled = x;
    end
end
xevaledout=xevaled;

xevaled = xevaled(:);
if model.SimpleLinearObjective
    f = model.f + model.c'*xevaled;
else
    f = model.f + (model.c'+xevaled'*model.Q)*xevaled;
    if isnan(f)
        f = inf;
    end
end
f=full(f);
df = [];
dx = [];
if nargout==1 || ~model.derivative_available
    return
elseif model.SimpleLinearObjective
    df = full(model.c(model.linearindicies));
elseif model.SimpleQuadraticObjective
    df = full(model.c(model.linearindicies) + 2*model.Q(model.linearindicies,model.linearindicies)*x);
elseif model.SimpleNonlinearObjective
    requested = model.c | any(model.Q,2);
    [i,j,k] = find((model.deppattern(find(requested),:)));
    requested(j) = 1;
    df = [];
    n = length(model.c);
    linearindicies = model.linearindicies;
    mtNonlinear = model.monomtable(model.nonlinearindicies,:);
    xevaled = zeros(1,n);
    xevaled(linearindicies) = x;
    X = repmat(xevaled,size(mtNonlinear,1),1);
    r = find(mtNonlinear);
    used = find(any(mtNonlinear,1));
    mtCompressed = mtNonlinear(:,used);
    rCompressed = find(mtCompressed);
    X = X(r);
    Xones = ones(size(mtNonlinear,1),size(mtCompressed,2));
    for i = 1:length(linearindicies)
        if requested(i)
            mt = mtNonlinear;
            oldpower = mtNonlinear(:,linearindicies(i));
            mt(:,linearindicies(i)) = mt(:,linearindicies(i))-1;
            Z = X.^mt(r);
            XX = Xones;          
            XX(rCompressed) = Z;
            xevaledNonLinear = prod(XX,2);
            xevaledNonLinear = xevaledNonLinear(:)'.*oldpower';xevaledNonLinear(isnan(xevaledNonLinear))=0;
            dx = zeros(1,n);
            dx(linearindicies(i)) = 1;
            dx(model.nonlinearindicies) = xevaledNonLinear;
            df = [df;model.c'*dx'];
        else
            df = [df;zeros(1,length(n))];
        end
    end
    df = real(df + 2*model.Q(model.linearindicies,model.linearindicies)*x);
    df = full(df);
elseif nargout > 1
    requested = model.c | any(model.Q,2);
    [i,j,k] = find((model.deppattern(find(requested),:)));
    requested(j) = 1;
    dx = apply_recursive_differentiation(model,xevaled,requested,model.frecursivederivativeprecompute); 
    df = model.c'*dx+2*xevaled'*model.Q*dx;
    df = full(df);
else
    df = [];
end






