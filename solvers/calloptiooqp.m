function output = calloptiooqp(interfacedata)

% Author Johan Löfberg 
% $Id: callopticlp.m,v 1.6 2005-05-07 13:53:20 joloef Exp $

% Standard input interface
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = full(interfacedata.c);
K       = interfacedata.K;
lb      = full(interfacedata.lb);
ub      = full(interfacedata.ub);
Q       = interfacedata.Q;

showprogress('Calling OOQP',options.showprogress);

blow = -inf(K.l,1);
if isempty(F_struc)  
    A = sparse([]);
    bupp = [];
    blow = [];
    Aeq = sparse([]);
    beq = [];
else
    Aeq = -F_struc(1:K.f,2:end);
    beq = full(F_struc(1:K.f,1));
    A =-F_struc(K.f+1:end,2:end);
    bupp = full(F_struc(K.f+1:end,1));
end

if isempty(lb)
    lb = -inf(length(c),1);
    ub = inf(length(c),1);
end

opts = options.ooqp;

H = 2*sparse(triu(Q));
A = A';
Aeq = Aeq';
if K.f == 0
    Aeq = [];
    beq = [];
end
if K.l == 0
    A = [];
    blow = [];
    bupp = [];
end
if options.savedebug
    save ooqpdebug c H lb ub Aeq beq A blow bupp opts
end

solvertime = clock; 
[x,fval,stat,iter] = ooqp(H, c, A, blow, bupp, Aeq,beq,lb,ub,opts);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

% No duals
D_struc = [];

%   Return Status:
%       0 - optimal
%       1 - not finished
%       2 - maximum iterations exceeded
%       3 - infeasible
%       4 - ooqp error 
switch stat
    case 0
        problem = 0;
    case 2
        problem = 3;
    case 3
        problem = 1;
    case 1
        problem = 11; 
    case {1,4}
        problem = 9;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,'OOQP');       

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.f = c;
    solverinput.lb = lb;
    solverinput.ub = ub;    
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.flag = flag;
    solveroutput.output=output;
else
    solveroutput = [];
end

% Standard interface 
output.Primal      = x(:);
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;