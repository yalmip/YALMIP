function output = kktqp(interfacedata)
%KKTQP Solver for indefinite QP problems using binary optmization
%
% Note that this solver can be implemented using very little code by using
% the high-level KKT operator. This version is only kept for compatibility
% reasons.
%
% See also kkt


% Author Johan Löfberg
% $Id: kktqp.m,v 1.2 2008-06-10 14:47:59 joloef Exp $

% ********************************
%% INITIALIZE DIAGNOSTICS IN YALMIP
% ********************************
K = interfacedata.K;
b = interfacedata.F_struc(1+K.f:end,1);
A = -interfacedata.F_struc(1+K.f:end,2:end);
c = interfacedata.c;

n = length(c); 

if ~all(isinf(interfacedata.lb))
    if ~isempty(interfacedata.lb)
        A = [A;-eye(n)];
        b = [b;-interfacedata.lb];
    end
end
if ~all(isinf(interfacedata.ub))
    if ~isempty(interfacedata.ub)
        A = [A;eye(n)];
        b = [b;interfacedata.ub];        
    end
end

x = sdpvar(length(c),1);
[dummy, L, U,sols,vals] = boundingbox(A*x <= b,[],[],interfacedata.f + interfacedata.c'*x + x'*interfacedata.Q'*x);

% Formulation here assumes maximization...
Q = -2*interfacedata.Q;
c = -interfacedata.c;

yalmip('setbounds',getvariables(x),L,U);

y = sdpvar(length(b),1);  % Duals 
dy = binvar(length(b),1); % indicater dual ==0 
ds = binvar(length(b),1); % indicater slack==0 
s = b-A*x;                % slack 

% Derive bounds on primal slack
[M,m] = derivebounds(s);

% Let us try to derive bounds on the dual variables
F = [A'*y == Q*x + c, s>=0, y>=0]; % KKT 
F = [F, s <= ds.*M];               % Big M, we know upper bound on s 
F = [F, dy+ds <= 1];               % Complementary slackness 
F = [F, 0 <= sum(dy) <= n];        % No need to force more than n active

% Find dis-joint constraints (silly way...)
for i = 1:length(b)
    j = findrows(abs(A),abs(A(i,:)));
    if length(j)>1
        S{i} = setdiff(j,i);
    else
        S{i} = [];
    end
end

if ~isempty(vals)
    F = [F, -0.5*(c'*x+b'*y) <= min([vals{:}])];
end

solvertime = clock;
My = derivedualbounds(F,y,b,S,n);
F = F + (y <= dy.*My);

obj = -0.5*(c'*x+b'*y); % ==cost in optimal points 
sol = solvesdp(F,obj,sdpsettings(interfacedata.options,'solver',interfacedata.options.kktqp.solver));

% **********************************
%% CREATE SOLUTION
% **********************************
output.problem = sol.problem;
output.Primal      = double(x);
output.Dual        = [];
output.Slack       = [];
output.infostr      = yalmiperror(output.problem,'KKTQP');
output.solverinput  = 0;
if interfacedata.options.savesolveroutput
    output.solveroutput.solved_nodes = solved_nodes;
    output.solveroutput.lower = lower;
    output.solveroutput.upper = upper;    
else
    output.solveroutput =[];
end
output.solvertime = etime(clock,solvertime);

function My = derivedualbounds(F,y,b,S,n)

[a1,a2,a3,model] = export(F,-y(1),sdpsettings('relax',1,'verbose',0),0,0,0);
for i = 1:length(b)
    if isempty(S{i})
        model_ = model;
        model_.c  = model_.c*0;
        model_.c(n+i)  = -1;
        sol = feval(model.solver.call,model_);
    else
        k = length(S{i});
        Aeq = sparse(1:k,n+S{i},1,k,n+3*length(b));
        beq = zeros(k,1);
        model_ = model;
        model_.F_struc = [beq Aeq;model.F_struc];
        model_.K.f = model.K.f+k;
        model_.c  = model_.c*0;
        model_.c(n+i)  = -1;
        sol = feval(model.solver.call,model_);
    end
    if sol.problem == 0
        My(i,1) =  sol.Primal(n+i); 
        model.ub(n+i) = My(i,1);
    else
        My(i,1) =  1e3;
    end
end
