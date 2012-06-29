function output = callbpmpd(interfacedata)

% Author Johan Löfberg 
% $Id: callbpmpd.m,v 1.4 2005-05-07 13:53:20 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;
 
if ~isempty(F_struc)
    A =-F_struc(1:end,2:end);
    b = F_struc(1:end,1);  
    e = [zeros(K.f,1);-ones(K.l,1)];
else
    A = [];
    b = [];
    e = [];
end

if isempty(ub) | all(isinf(lb))
    llist = [];
    lval = [];
else
    llist = find(~isinf(lb));
    lval = lb(llist);       
end

if isempty(ub) | all(isinf(ub))
    ulist = [];
    uval = [];
else
    ulist = find(~isinf(ub));
    uval = ub(ulist);   
end

Q = 2*Q;
b = full(b);
c = full(c);

if nnz(Q)==0
    Q = [];
end

opts = options.bpmpd.opts;
if options.savedebug
    ops = options.quadprog;
    save bpmpdebug Q A b c e llist lval ulist uval opts
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = clock; 
[x,y,s,w,how] = bp(Q, A, b, c, e,llist,lval,ulist,uval,opts);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

problem = 0;

% Internal format for duals
D_struc = -y;

switch how
    case 'optimal solution'
        problem = 0;
    case 'infeasible primal'
        problem = 1;
    case 'infeasible dual'
        problem = 2;
    otherwise
        problem = 11;
end

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.c = c;
    solverinput.Q = Q;
    solverinput.e = e;
    solverinput.llist = llist;
    solverinput.lval = lval;
    solverinput.ulist = ulist;
    solverinput.uval = uval;    
    solverinput.bpoptss = bpopt;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.y = y;
    solveroutput.s = s;
    solveroutput.w = w;
    solveroutput.how=how;  
else
    solveroutput = [];
end



% Standard interface 
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;