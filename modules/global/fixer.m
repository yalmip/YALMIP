function [upper1,x_min] = fixer(p,output);

x = output.Primal;

integer = uniquestripped([p.binary_variables p.integer_variables]);
n = length(integer);
Fnew = [round(x(integer)) -sparse(1:n,integer,1,n,length(p.c))];
p.F_struc = [Fnew;p.F_struc];
p.K.f = p.K.f + n;

% close = find(abs(round(x(integer)) - x(integer))<1e-2);
% n = length(integer(close));
% Fnew = [round(x(integer(close))) -sparse(1:n,integer(close),1,n,length(p.c))];
% p.F_struc = [Fnew;p.F_struc];
% p.K.f = p.K.f + n;
% % 

% Add the bounds
upper = find(~isinf(p.ub));
if ~isempty(upper)
    Fnew = [p.ub(upper) -sparse(1:length(upper),upper,1,length(upper),length(p.c))];
    p.F_struc = [p.F_struc(1:p.K.f,:);Fnew;p.F_struc(p.K.f + 1:end,:)];
    p.K.l = p.K.l + length(upper);
end
lower = find(~isinf(p.lb));
if ~isempty(lower)
    Fnew = [-p.lb(lower) sparse(1:length(lower),lower,1,length(lower),length(p.c))];
    p.F_struc = [p.F_struc(1:p.K.f,:);Fnew;p.F_struc(p.K.f + 1:end,:)];
    p.K.l = p.K.l + length(lower);
end
p.ub = p.ub + inf;
p.lb = p.lb - inf;
try
    output = bnb_solvelower(p.solver.lower.call,p,inf,nan);
    x  = setnonlinearvariables(p,output.Primal);
    if output.problem == 0
        x_min  = setnonlinearvariables(p,output.Primal);        
        %x_min = output.Primal;
         upper1 = computecost(p.f,p.c,p.Q,x_min,p);%upper1 = x_min'*p.Q*x_min + p.c'*x_min + p.f;
    else
        x_min = [];
        upper1 = inf;
    end
catch
    x_min = [];
    upper1 = inf;
end

