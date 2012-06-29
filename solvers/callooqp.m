function output = callquadprog(interfacedata)

% Author Johan Löfberg 
% $Id: callooqp.m,v 1.3 2006-11-29 16:45:58 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;

switch options.verbose
case 0
    options.quadprog.Display = 'off';
case 1
    options.quadprog.Display = 'final';
otherwise
    options.quadprog.Display = 'iter';
end
    
showprogress('Calling QUADPROG',options.showprogress);

if ~isempty(F_struc)
    A = -F_struc(1:1:K.f,2:end);
    dA = full(F_struc(1:1:K.f,1));        
    C =-F_struc(K.f+1:end,2:end);
    cupp = full(F_struc(K.f+1:end,1));   
    clow = full(-ones(size(C,1),1)*1e4);
else
    A  = [];
    dA = [];
    C  = [];
    cupp = [];
    clow = [];
end

n = length(c);
if isempty(lb)
    xlow = -ones(n,1)*1e9;
    xupp = ones(n,1)*1e9;
else
    lb(isinf(lb)) = -1e9;
    ub(isinf(ub)) = 1e9;
    xlow = lb;
    xupp = ub;
end

if options.verbose==0
    doprint = 'no';
else
    doprint = 'yes';
end
solvertime = clock; 
[status, x, gamma, phi, y, z, lambda, pi] = ooqp( c, 2*Q, xlow, xupp, A, dA, C, clow, cupp, doprint);
solvertime = etime(clock,solvertime);

problem = 0;

% Internal format for duals
D_struc = -[y;z];

% Check, currently not exhaustive...
switch status
    case 0
        problem = 0;
    case 3
        problem = 1;
    otherwise
        problem = 13;
end
infostr = yalmiperror(problem,'OOQP');       

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.Aeq = Aq;
    solverinput.beq = beq;
    solverinput.c = c;
    solverinput.H = Q;
    solverinput.options = options.quadprog;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.flag = flag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;  
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