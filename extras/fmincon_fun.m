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
    [the_linearindex,the_nonlinearindex,the_powers] = find(mtNonlinear');
    if isequal(the_nonlinearindex(:)',1:size(mtNonlinear,1))
        trivialUnivariateCase = 1;
    else
        trivialUnivariateCase = 0;
    end
    xpower_all = x(the_linearindex).^the_powers;
    can_nan_occur = any(model.variabletype > 3);
    for i = 1:length(linearindicies)
        if requested(i)            
            oldpower = mtNonlinear(:,linearindicies(i));
                        
            xpower = xpower_all;
            currentdifferentiation = find((the_linearindex == i));
            xpower(currentdifferentiation) = x(the_linearindex(currentdifferentiation)).^(the_powers(currentdifferentiation)-1);
            if trivialUnivariateCase
                product = xpower.*(oldpower);
            else
                product = accumarray(the_nonlinearindex,xpower,[],@prod).*(oldpower);
            end
            xevaledNonLinear = product(:)';                    
            if can_nan_occur
                xevaledNonLinear(isnan(xevaledNonLinear)) = 0;
            end                      
            dx = zeros(1,n);
            dx(linearindicies(i)) = 1;
            dx(model.nonlinearindicies) = xevaledNonLinear;
            c_x = model.c'*dx';            
            df = [df;c_x];
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






