function [dualUpper,L,U] = derivedualBoundsParameterFree(H,c,A,b,E,f,ops,parametricDomain)

if isempty(A)
    dualUpper = [];
    L = [];
    U = [];
    return
end

n = length(c);
m = length(b);
me = length(f);

x = sdpvar(n,1);
Lambda = sdpvar(m,1);
mu = sdpvar(me,1);

% Homogenization
alpha = sdpvar(1);

F = [];

% Start by computing primal lower bounds
if nargin < 7
    ops = sdpsettings('verbose',0);
end
ops2 = ops;
ops2.verbose = max(0,ops.verbose-1);;
all_bounded = 1;
if ops.verbose
    disp(['*Computing ' num2str(length(x)) ' primal bounds (required for dual bounds)']);
end

z = recover(unique([depends(c);depends(b)]));
xz = [x;z];
nz = length(z);
nTOT = n + length(z);
rhs = 0;
if ~isempty(b)
    rhs =  rhs + A'*Lambda;
end
if ~isempty(f)
    rhs =  rhs + E'*mu;
end

[c0,C] = deParameterize(c,z);
[b0,B] = deParameterize(b,z);
delta = binvar(length(b),1);
UpperBound = .1;
parametricDomainH = homogenize(sdpvar(parametricDomain),alpha) >= 0;

w1 = binvar(length(Lambda),1);
w2 = binvar(length(Lambda),1);
eta1 = sdpvar(length(Lambda),1);
eta2 = sdpvar(length(Lambda),1);
v = sdpvar(size(A,2),1);
Model = [parametricDomainH,H*x + alpha*c0 + C*z + rhs==0,
                           delta >= Lambda >= 0,
                         1-delta >= b0*alpha+B*z-A*x>=0, 
                         sum(Lambda) >= alpha*UpperBound,
                       %  sum(delta)<=1,
                         ];
Model = [Model, Lambda + A*v + eta2-eta1 == 0,
                 0 <= eta1 <= w1, 0 <= delta - Lambda <= 1-w1,
                 0 <= eta2 <= w2, 0 <= Lambda <= 1-w2];
                               
solvesdp(Model,-alpha)
while double(alpha) >= 1e-5
    UpperBound = UpperBound*1.1;
    Model = [parametricDomainH,H*x + alpha*c0 + C*z + rhs==0,
             delta >= Lambda >= 0,
             1-delta >= b0*alpha+B*z-A*x>=0,
             sum(Lambda) >= alpha*UpperBound];
          Model = [Model, Lambda + A*v + eta2-eta1 == 0,
                  0 <= eta1 <= w1, 0 <= delta - Lambda <= 1-w1,
                  0 <= eta2 <= w2, 0 <= Lambda <= 1-w2];
    solvesdp(Model,-alpha);
end  

dualUpper = min(U,UpperBound);    



function [c0,C] = deParameterize(c,z);
n = length(c);
c0C = full(getbase(c));
c0 = c0C(:,1);
C = c0C(:,2:end);
if nnz(C)==0
    C = spalloc(n,length(z),0);
else
    if length(z)>0 & size(C,2) < length(z)
        C = [];
        for i = 1:length(z)
            C = [C getbasematrix(c,getvariables(z(i)))];
        end  
    end
end

if isempty(C)
    C = zeros(length(C),length(z));
end