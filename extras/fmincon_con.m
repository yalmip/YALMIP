function [g,geq,dg,dgeq,xevaled] = fmincon_con(x,model,xevaled)

global latest_xevaled
global latest_x_xevaled
global sdpLayer

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
        [d,v] = eig(X);        
        v = sdpLayer.f(diag(v));
        % These will reordered later
        g = [g;-v(1:min(n,sdpLayer.n))];
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
    dg = computeSOCPDeriv(model,xevaled);
    % Linear SDP
    dg_SDP = computeSDPDeriv(model,xevaled,1);
    dg = [dg;dg_SDP];     
    dg = dg(:,model.linearindicies);      
    g = reorderEigenvalueG(g,model); 
    
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
        dg = [dg;computeSOCPDeriv(model,xevaled,newdxx)];
    end          
    if any(model.K.s)
        dg_SDP = computeSDPDeriv(model,xevaled,newdxx);               
        dg = [dg;dg_SDP]; 
        g = reorderEigenvalueG(g,model);       
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
        dg = [dg;computeSOCPDeriv(model,xevaled,dx)];
    end    
    if any(model.K.s)
        dg_SDP = computeSDPDeriv(model,xevaled,dx);                       
        dg = [dg;dg_SDP];
        g = reorderEigenvalueG(g,model);       
    end
end

% For performance reasons, we work with transpose
if model.nonlinearequalities
    dgeq = dgeq';
end
if model.nonlinearinequalities | any(model.K.q) | any(model.K.s)
    dg = dg';   
end
%  x = sdpvar(2,1);
%  plot(g + dg'*(x-xevaled(1:2))<=0,x,[],[],sdpsettings('plot.shade',.1))
% plot(xevaled(1),xevaled(2),'+b')
% drawnow
% axis([-2 2 -2 2])
% 1

function conederiv = computeSOCPDeriv(model,z,dzdx)
conederiv = [];
z = sparse(z(:));
if any(model.K.q)
    top = 1 + model.K.f + model.K.l;  
    for i = 1:length(model.K.q)
        d = model.F_struc(top,1);
        b = model.F_struc(top+1:top+model.K.q(i)-1,1);
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
        else
            % -c'*x - d + ||Ax+b||>=0
            e = A*z + b;
            
            % smoothed = sqrt(10^-10 + e'*e);
            % aux = z'*(A'*A-c*c')*dzdx+(b'*A-d*c')*dzdx;
            % conederiv = [conederiv ((-dzdx'*c+aux'/smoothed)')'];

            smoothed = sqrt(10^-10 + e'*e);
            temp = (-c'+(A'*(b + A*z))'/smoothed);
            conederiv = [conederiv (temp*dzdx)'];
        end
        top = top + model.K.q(i);
    end
end
conederiv = conederiv';

function [dg,reordering] = computeSDPDeriv(model,xevaled,newdxx)
top = 1 + sum(model.K.q);
newF = [];newcuts = 1;
B=model.F_struc(:,2:end)*newdxx;
global sdpLayer
reordering = [];

for i = 1:length(model.K.s)
    n = model.K.s(i);
    X = model.F_struc(top:top+n^2-1,:)*[1;xevaled];
    X = full(reshape(X,n,n));
    [d,v] = eig(X);    
    newSDPblock = [];
    for m = 1:min(sdpLayer.n,model.K.s(i))
        newrow = [];
        for j = 1:size(B,2)
            newrow = [newrow d(:,m)'*reshape(B(top:top+n^2-1,j),n,n)*d(:,m)];
        end
       % newrow = reshape(d(:,m)*d(:,m)',[],1)'*B(top:top+n^2-1,:);
        newrow = newrow*sdpLayer.df(v(m,m));
        newSDPblock = [newSDPblock;newrow];  
        newcuts = newcuts + 1;
    end     
    if ~isempty(sdpLayer.oldGradient{i}) 
        [reordering] = matchGradientRows(newSDPblock,sdpLayer.oldGradient{i});
        newSDPblock = newSDPblock(reordering,:);                        
    end
    sdpLayer.oldGradient{i} = newSDPblock;
    sdpLayer.reordering{i} = reordering;
    newF = [newF;newSDPblock];
    top = top + n^2;
end
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

function r = matchGradientRows(X,Y)

C = distancematrix(X,Y);
[r,m,u] = matchpairs(C,1e12);
r = r(:,1);

function C = distancematrix(X,Y)

C = zeros(size(X,1));
for i = 1:size(X,1)
    for j = 1:size(X,1)
        C(i,j) = norm(X(i,:)-Y(j,:),1);               
    end
end

function  g = reorderEigenvalueG(g,model);
global sdpLayer
if ~isempty(sdpLayer.reordering{1})
    g_nonsdp = g(1:end-sum(min(sdpLayer.n,model.K.s)));
    g_sdp = g(end-sum(min(sdpLayer.n,model.K.s))+1:end);
    reordering = [];
    top = 0;
    for i = 1:length(model.K.s)
        reordering = [reordering;top+sdpLayer.reordering{i}];
        top = top+ min(sdpLayer.n,model.K.s(i));
    end
    g_sdp = g_sdp(reordering);
    g = [g_nonsdp;g_sdp];
end
