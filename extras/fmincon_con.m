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

if nnz(model.K.q) > 0
    top = 1;
    for i = 1:length(model.K.q)
        z = model.F_struc(top:top+model.K.q(i)-1,:)*[1;xevaled];
        g = [g;-(z(1)^2 - z(2:end)'*z(2:end))];
        top = top + model.K.q(i);
    end
end

if model.nonlinearequalities
    geq = full(model.Anonlineq*xevaled(:)-model.bnonlineq);
end

dgAll_test = [];

if nargout == 2 || ~model.derivative_available
    return
elseif ~isempty(dgAll_test) & isempty(model.evalMap) 
    dgAll = dgAll_test;
elseif isempty(model.evalMap) & (model.nonlinearinequalities==0) & (model.nonlinearequalities==0) & (model.nonlinearcones==0) & any(model.K.q)
    dg = computeConeDeriv(model,xevaled);
elseif isempty(model.evalMap) & (model.nonlinearinequalities | model.nonlinearequalities | model.nonlinearcones) 
    n = length(model.c);
    linearindicies = model.linearindicies;    
 %   xevaled = zeros(1,n);
 %   xevaled(linearindicies) = x;
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
    %newdxx = newdxx';
    
    if ~isempty(model.Anonlineq)      
        dgeq = model.Anonlineq*newdxx; 
    end
    if ~isempty(model.Anonlinineq)
        dg = model.Anonlinineq*newdxx;        
    end
    
    if nnz(model.K.q)>0
        dg = [dg;computeConeDeriv(model,xevaled,newdxx);];
    end    
else    
    requested = model.fastdiff.requested;
    dx = apply_recursive_differentiation(model,xevaled,requested,model.Crecursivederivativeprecompute);    
    conederiv = computeConeDeriv(model,xevaled,dx);   
    if ~isempty(model.Anonlineq)
        dgeq = [model.Anonlineq*dx];  
    end
    if ~isempty(model.Anonlinineq)
        dg = [model.Anonlinineq*dx];
    end    
    if ~isempty(conederiv)
        dg = [dg;conederiv];
    end
end

if model.nonlinearequalities
    dgeq = dgeq';
end
if model.nonlinearinequalities | any(model.K.q)
    dg = dg';
end
    

function conederiv = computeConeDeriv(model,z,dzdx)
conederiv = [];
z = z(:);
if any(model.K.q)
    top = 1 + model.K.f + model.K.l;  
    for i = 1:length(model.K.q)
        d = model.F_struc(top,1);
        c = model.F_struc(top,2:end)';
        b = model.F_struc(top+1:top+model.K.q(i)-1,1);
        A = model.F_struc(top+1:top+model.K.q(i)-1,2:end);
        
        if nargin == 2
            % No inner derivative
            
            conederiv = [conederiv;(2*A(:,model.linearindicies)'*(A(:,model.linearindicies)*z(model.linearindicies)+b)-2*c(model.linearindicies)*(c(model.linearindicies)'*z(model.linearindicies)+d))'];
        else
            % inner derivative
            aux = 2*z'*(A'*A-c*c')*dzdx+2*(b'*A-d*c')*dzdx;
            conederiv = [conederiv;aux];
        end
        top = top + model.K.q(i);
    end
end
