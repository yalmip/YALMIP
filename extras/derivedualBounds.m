function [dualUpper,L,U] = derivedualBounds(H,c,A,b,E,f,ops,parametricDomain)

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

F = [];
if ~isempty(b)
    F = [F, A*x <= b];
end
if ~isempty(f)
    F = [F, E*x == f];
end

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

z = recover(unique([depends(c) depends(b)]));
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

[dummy,L,U] = boundingbox([F,parametricDomain,H*x + c + rhs==0,Lambda>=0],ops2,xz);

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

b0B = full(getbase(b));
b0 = b0B(:,1);
B = b0B(:,2:end);
if nnz(B)==0
    B = spalloc(n,length(z),0);
else
    if length(z)>0 & size(B,2) < length(z)
        B = [];
        for i = 1:length(z)
            B = [B getbasematrix(b,getvariables(z(i)))];
        end
    end
end
    
if isa(b,'sdpvar')
    if ops.verbose
        disp(['*Computing ' num2str(length(b)) ' bounds on parameterized RHS (required for dual bounds)']);
    end  
    [dummy,bmin,bmax] = boundingbox([F,parametricDomain],ops2,b);
end

% Lift
Hlift = [H C/2;C'/2 zeros(size(C,2))];
if ops.kkt.dualpresolve.lplift
    XZ = sdpvar(nTOT);
end
constraints = [Lambda >= 0,parametricDomain,F];
constraints = [constraints,H*x + c + rhs==0];
finiteL = find(~isinf(L));
if ~isempty(finiteL)
    constraints = [constraints,L(finiteL)<=xz(finiteL)];
end
finiteU = find(~isinf(U));
if ~isempty(finiteU)
    constraints = [constraints,xz(finiteU) <= U(finiteU)];
end

if ~isa(f,'double')
    disp('Parameterized RHS in equality not yet supported in lifting based bounds')
    ops.kkt.dualpresolve.lplift = 0;
end

if ops.kkt.dualpresolve.lplift
    if isempty(f)
        temp = 0;
    else
        temp = f'*mu;
    end
           
    if nnz(Hlift)>0
        for i = 1:nTOT
            for j = i:nTOT 
                if (L(i)==U(i) & L(j)==U(j)) & ~isinf(L(j))               
                   XZ(i,j) = L(i)*L(j);
                   XZ(j,i) = L(i)*L(j);
                   XZ(i,i) = L(i)^2;
                   XZ(j,j) = L(j)^2;
                elseif L(i)==U(i) & ~isinf(L(i))
                    XZ(i,j) = L(i)*xz(j);
                    XZ(j,i) = L(i)*xz(j);
                    XZ(i,i) = L(i)^2;
                elseif L(j)==U(j) & ~isinf(L(j))
                    XZ(i,j) = L(j)*xz(i);
                    XZ(j,i) = L(j)*xz(i);            
                    XZ(j,j) = L(j)^2;
                else
                    % (xi-Li)(xj-Lj)>0                    
                    if ~isinf(L(i)) & ~isinf(L(j))
                        constraints = [constraints, XZ(i,j)-L(i)*xz(j)-L(j)*xz(i)+L(i)*L(j)>=0];
                    end
                    % (xi-Li)(Uj-xj)>0
                    if ~isinf(L(i)) & ~isinf(U(j))
                        constraints = [constraints, xz(i)*U(j)+L(i)*xz(j)-L(i)*U(j)-XZ(i,j)>=0];
                    end
                    % (xj-Lj)(Ui-xi)>0
                    if ~isinf(U(i)) & ~isinf(L(j))
                        constraints = [constraints, xz(j)*U(i)+L(j)*xz(i)-L(j)*U(i)-XZ(j,i)>=0];
                    end
                    % (Uj-xj)(Ui-xi)>0
                    if ~isinf(U(i)) & ~isinf(U(j))
                        constraints = [constraints, U(j)*U(i)-U(j)*xz(i)-xz(j)*U(i)+XZ(j,i)>=0];
                    end
                end
            end
        end
    end
    
    if isa(b,'sdpvar')
        if ~any(isinf(bmax))
            constraints = [constraints,-(trace(Hlift*XZ) + c0'*x+temp) <= bmax'*Lambda];
        end
        if ~any(isinf(bmin))
            constraints = [constraints,-(trace(Hlift*XZ) + c0'*x+temp) >= bmin'*Lambda];
        end
    else
        constraints = [constraints,trace(Hlift*XZ) + c0'*x + temp + b'*Lambda==0];
    end
    
end

if ops.verbose
    disp(['*Computing ' num2str(length(Lambda)) ' dual bounds']);
end
ops2.verbose = 0;

dualUpper = inf(length(Lambda),1);

% always exploit disjoint constraints
disjoints = zeros(length(b));
s = sdpvar(length(b),1);

% Create an optimizer object which takes no input and returns optimal
% Lambda. We will use some hacks to update this model and resolving
% with different objectives and constraints. We cannot use multiple
% objectives approach since we have different constraints, and we
% cannot use optimizer object directly since constraints change
p = optimizer([constraints,s == b-A*x],-sum(Lambda),ops2,[],Lambda);

Model = p.model;
LambdaIndicies = find(Model.c);
Model.c = Model.c*0;
sIndicies = [ find(Model.used_variables==getvariables(s(1))): find(Model.used_variables==getvariables(s(end)))];

for dualPass = 1:ops.kkt.dualpresolve.passes    
    for i = 1:m
        if dualUpper(i)>0           
            iModel = Model;
            % Maximize Lambda(i)...
            iModel.c = Model.c*0;
            iModel.c(LambdaIndicies(i))=-1;
            
            % s(i) == 0
            iModel.lb(sIndicies(i)) = 0;
            iModel.ub(sIndicies(i)) = 0;
                        
            % Add known dual bounds
            r = find(~isinf(dualUpper));

            iModel.ub(LambdaIndicies(r)) = dualUpper(r);
            
            % Find disjoint constraints
            a = A(i,:);
            a = setdiff(findrows(A,-a),i);
            iModel.ub(LambdaIndicies(a)) = 0;
            
            % Solve this instance
            p = updatemodel(p,iModel);
            [Lambda,isol] = p{[]};
            sol.problem = isol;
            
            if sol.problem == 12
                iModel.c = iModel.c*0;
                p = updatemodel(p,iModel);
                [Lambda,isol] = p{[]};
                sol.problem = isol;
                if sol.problem == 1
                    dualUpper(i,1) = 0;
                elseif sol.problem == 0
                    sol.problem = 2;
                end
            end
             if sol.problem == 0
                dualUpper(i,1) = double(Lambda(i));
            elseif sol.problem == 2
                dualUpper(i,1) = inf;
            end
        end
    end
end

if ops.verbose
    if all(isinf(dualUpper))
        disp('*Warning: All dual upper bounds are infinite');
    elseif any(isinf(dualUpper))
        disp('*Warning: Some of the dual upper bounds are infinite');
    else
        disp('*Success: All dual upper bounds are finite');
    end
end




