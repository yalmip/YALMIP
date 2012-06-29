function [upper,x_min] = feaspump(p,x)

% Assume we fail
upper = inf;
x_min = x;

% These should be integer
intvars = [p.integer_variables(:);p.binary_variables(:)];
n = length(p.c);

% Create rounded goal variable
xround = x;
xround(intvars) =  round(xround(intvars));

p.c = [0*p.c;ones(n,1)];
p.Q = blkdiag(0*p.Q,spalloc(n,n,0));

f1 = [p.F_struc(1:p.K.f,:) zeros(p.K.f,n)];
f2 = [p.F_struc(1+p.K.f:end,:) zeros(size(p.F_struc,1)-p.K.f,n)];
p.F_struc = [f1;xround -eye(n) eye(n);-xround eye(n) eye(n);f2];
p.K.l = p.K.l + 2*n;
p.lb = [p.lb;-inf(n,1)];
p.ub = [p.ub;inf(n,1)];

% Append model with absolute value terms
c = p.c;
Q = p.Q;
p.c = [0*p.c;ones(n,1)];
p.Q = blkdiag(0*p.Q,spalloc(n,n,0));
p.F_struc = [p.F_struc(1:p.K.f,:) zeros(p.K.f,n);
    xround -eye(n) eye(n);
    -xround eye(n) eye(n);
    p.F_struc zeros(size(p.F_struc,1)-p.K.f,n)];
p.monomtable = blkdiag(p.monomtable,eye(n));
p.variabletype = [p.variabletype zeros(1,n)];
p.K.l = p.K.l + 2*n;
p.lb = [p.lb;-inf(n,1)];
p.ub = [p.ub;inf(n,1)];

iter = 1;
xold = x;
while iter < 10 & isinf(upper)

    % Append model with absolute value terms
    ii = p.integer_variables;
    jj = p.binary_variables;
    p.integer_variables = [];
    p.binary_variables = [];
    p.F_struc(1+p.K.f:p.K.f+2*n,1) = [xround;-xround];
       
    pumpsol = feval(p.solver.lower.call,p);
    p.integer_variables = ii;
    p.binary_variables = jj;
    
    % Fix the varying data
    p.F_struc(1+p.K.f:1+p.K.f+2*n-1,1) = [xround;-xround];
    
    % Solve pumping problem
    pumpsol = feval(p.solver.lower.call,p);
    
    % New solution (extract original variables)
    xnew = pumpsol.Primal(1:n);
    
    if all(abs(xnew(intvars)-round(xnew(intvars))) < 1e-5)
        % Feasible
        upper = c'*xnew + xnew'*Q*xnew;
        x_min = xnew(1:n);
    else
        % Stall?
        if norm(xold-xnew) < 1e-3
            % Yep, so flip the goal on some integers
            [ii,jj] = sort(abs(xnew(intvars)-xround(intvars)));
            flips = jj(end-(1+ceil(5*rand(1))):end);
          %  for i  =1:length(jj)
                xround(intvars(flips)) = (xround(intvars(flips))-1).^2;
          %  end
        else
            xround = xnew;xold = xnew;
            xround(intvars) =  round(xround(intvars));
        end
    end
    iter = iter + 1;
end
