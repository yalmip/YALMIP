function [g,geq,dg,dgeq,xevaled] = fmincon_con(x,model,xevaled)

global latest_xevaled
global latest_x_xevaled

% Early bail for linear problems
g    = [];
geq  = [];
dg   = [];
dgeq = [];

% Nothing nonlinear, and no cones
if model.linearconstraints && ~any(model.K.q) && ~any(model.K.s)
    xevaled = [];
    return
end

% Compute all nonlinear variables
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

% All simple scalar inequalities
if model.nonlinearinequalities
    g = full(model.Anonlinineq*xevaled(:)-model.bnonlinineq);    
end

% Append SOCP cones
if any(model.K.q)
    top = 1;
    z0 = model.F_struc*[1;xevaled];
    for i = 1:length(model.K.q)
        z = z0(top:top+model.K.q(i)-1);
        g = [g;-(z(1) - sqrt(z(2:end)'*z(2:end)))];
        top = top + model.K.q(i);
    end
end

% Append SDP cones (eigenvalues)
if any(model.K.s)
    top = 1 + sum(model.K.q);
    for i = 1:length(model.K.s)
        n = model.K.s(i);
        X = model.F_struc(top:top+n^2-1,:)*[1;xevaled];
        X = full(reshape(X,n,n));
        d = eig(X);
        %d(d>1) = 0;
        g = [g;-d];
        top = top + n^2;
    end
end

% Create all nonlinear equalities
if model.nonlinearequalities
    geq = full(model.Anonlineq*xevaled(:)-model.bnonlineq);
end

% Now Jacobians...
if nargout == 2 || ~model.derivative_available
    return
elseif isempty(model.evalMap) & (model.nonlinearinequalities==0) & (model.nonlinearequalities==0) & (model.nonlinearcones==0) & (any(model.K.q) || any(model.K.s))
    % Linear SOCP
    dg = computeConeDeriv(model,xevaled);
    % Linear SDP
    dg = [dg;computeSDPDeriv(model,xevaled,1)];       
elseif isempty(model.evalMap) & (model.nonlinearinequalities | model.nonlinearequalities | model.nonlinearcones) 
    % Nonlinear terms in constraints so jacobian stuff needed
    % however, there are only monomials, which is exploited
    newdxx = computeMonomialVariableJacobian(model,x);
        
    if ~isempty(model.Anonlineq)      
        dgeq = model.Anonlineq*newdxx; 
    end
    if ~isempty(model.Anonlinineq)
        dg = model.Anonlinineq*newdxx;                
    end    
    if any(model.K.q)
        dg = [dg;computeConeDeriv(model,xevaled,newdxx)];
    end        
    if any(model.K.s)
        dg = [dg;computeSDPDeriv(model,xevaled,newdxx)];       
    end
else  
    % Completely general case with monomials and functions
    requested = model.fastdiff.requested;
    dx = apply_recursive_differentiation(model,xevaled,requested,model.Crecursivederivativeprecompute);        
    if ~isempty(model.Anonlineq)
        dgeq = [model.Anonlineq*dx];  
    end
    if ~isempty(model.Anonlinineq)
        dg = [model.Anonlinineq*dx];
    end        
    if any(model.K.q)
        dg = [dg;computeConeDeriv(model,xevaled,dx)];
    end    
    if any(model.K.s)
        dg = [dg;computeSDPDeriv(model,xevaled,dx)];
    end
end

% For performance reasons, we work with transpose
if model.nonlinearequalities
    dgeq = dgeq';
end
if model.nonlinearinequalities | any(model.K.q) | any(model.K.s)
    dg = dg';   
end
    

function conederiv = computeConeDeriv(model,z,dzdx)
conederiv = [];
z = sparse(z(:));
if any(model.K.q)
    top = 1 + model.K.f + model.K.l;  
    for i = 1:length(model.K.q)
        d = model.F_struc(top,1);
       % c = model.F_struc(top,2:end)';
        b = model.F_struc(top+1:top+model.K.q(i)-1,1);
       % A = model.F_struc(top+1:top+model.K.q(i)-1,2:end);
        cA = model.F_struc(top:top+model.K.q(i)-1,2:end);
        c = cA(1,:)';
        A = cA(2:end,:);
        
        if nargin == 2
            % No inner derivative
            if length(model.linearindicies)==length(model.c)
                e = A*z + b;
                smoothed = sqrt(10^-10 + e'*e);
                temp = -c+(A'*(b + A*z))/smoothed;
                conederiv = [conederiv temp];
            else
                A = A(:,model.linearindicies);
                c = c(model.linearindicies);
                % -c'*x - d + ||Ax+b||>=0
                e = A*z(model.linearindicies) + b;
                smoothed = sqrt(10^-10 + e'*e);
                temp = (-c'+(A'*(b + A*z(model.linearindicies)))'/smoothed);
                conederiv = [conederiv temp'];
            end
           % conederiv = [conederiv;(2*A(:,model.linearindicies)'*(A(:,model.linearindicies)*z(model.linearindicies)+b)-2*c(model.linearindicies)*(c(model.linearindicies)'*z(model.linearindicies)+d))'];
        else                       
            % -c'*x - d + ||Ax+b||>=0
            e = A*z + b;
            smoothed = sqrt(10^-10 + e'*e);
            % conederiv = [conederiv;-c'+(dzdx'*A'*b + dzdx'*A'*A*z(model.linearindicies))'/smoothed]; 
            aux = z'*(A'*A-c*c')*dzdx+(b'*A-d*c')*dzdx;
            conederiv = [conederiv ((-dzdx'*c+aux'/smoothed)')'];             
            % inner derivative
            % aux = 2*z'*(A'*A-c*c')*dzdx+2*(b'*A-d*c')*dzdx;
            % conederiv = [conederiv;aux];                                    
        end
        top = top + model.K.q(i);
    end
end
conederiv = conederiv';

function dg = computeSDPDeriv(model,xevaled,newdxx)
top = 1 + sum(model.K.q);
newF = [];newcuts = 1;
B=model.F_struc(:,2:end)*newdxx;
persistent xold

for i = 1:length(model.K.s)
    n = model.K.s(i);
    X = model.F_struc(top:top+n^2-1,:)*[1;xevaled];
    X = full(reshape(X,n,n));
    [d,v] = eig(X);  
    plot(xevaled(1),xevaled(2),'*');
    if 0%~isempty(xold)  && any(xevaled)              
        if 1
            sdpvar x1 x2
            s1 = sdisplay(d(:,1)'*reshape(full(model.F_struc)*[1;x1;x2],n,n)*d(:,1));
            try
                l=ezplot(s1{1});set(l,'LineColor','red');hold on
            catch
            end
            try
                s1 = sdisplay(d(:,2)'*reshape(full(model.F_struc)*[1;x1;x2;0],n,n)*d(:,2));
                ezplot(s1{1})
                plot(x(1),x(2),'k*');drawnow
            catch
            end
        end
    end
    xold = xevaled;
    X = model.F_struc(top:top+n^2-1,:)*[1;xevaled];
    X = full(reshape(X,n,n));
    [d,v] = eig(X);      
    for m = 1:n
        % sum(kron([1;1],v).*kron(v,[1;1]).*q)??
        %newF = [newF;reshape(d(:,m)*d(:,m)',[],1)'*model.F_struc(top:top+n^2-1,:)];             
        %newF = [newF;(v(m,m)<=1)*reshape(d(:,m)*d(:,m)',[],1)'*B(top:top+n^2-1,:)];        
         newF = [newF;reshape(d(:,m)*d(:,m)',[],1)'*B(top:top+n^2-1,:)];        
       % q1 = kron(ones(n,1),d(:,m));
       % q2 = kron(d(:,m),ones(n,1));
       % newF = [newF;sum((q1.*q2).*model.F_struc(top:top+n^2-1,:))];
        newcuts = newcuts + 1;
    end
    top = top + n^2;
end
%dg = -newF(:,2:end)*newdxx;
dg = -newF;

function newdxx = computeMonomialVariableJacobian(model,x)
n = length(model.c);
linearindicies = model.linearindicies;

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