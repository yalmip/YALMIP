function output = calllinprog(interfacedata)

% Author Johan Löfberg 
% $Id: callbintprog.m,v 1.5 2007-02-08 13:51:23 joloef Exp $

% Standard input interface
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;
x0      = interfacedata.x0;
binary_variables = interfacedata.binary_variables;
integer_variables = interfacedata.integer_variables;

showprogress('Calling BINTPROG',options.showprogress);

if isempty(F_struc)
    Aeq = [];
    beq = [];
    A = [];
    b = [];
else
    Aeq = -F_struc(1:1:K.f,2:end);
    beq = F_struc(1:1:K.f,1);        
    A =-F_struc(K.f+1:end,2:end);
    b = F_struc(K.f+1:end,1);   
end
solvertime = clock; 

% Any continuous varialbles
if length(union(binary_variables,integer_variables)) < length(c)
    % Standard interface
    output.Primal      = [];
    output.Dual        = [];
    output.Slack       = [];
    output.problem     = -4;
    output.infostr     = yalmiperror(-4);
    output.solverinput = [];
    output.solveroutput= [];
    output.solvertime  = 0;
    return
end

if any(integer_variables)
    [lb,ub,used_rows] = findulb(F_struc,K);
    % Not explicitely bounded not allowed
    if any(any(isinf([lb(integer_variables) ub(integer_variables)])))
        % Standard interface
        if options.verbose
            disp('Error when preparing call to BINTPROG.')
            disp('There are unbounded integer varibles in your model!')
            disp('BINTPROG can only handle binary variables, and when')
            disp('YALMIP converts integer variables to binary, they')
            disp('must have explicit bounds in the model');
            disp('Solution: use some other MILP solver (such as BNB)')
        end
        output.Primal      = [];
        output.Dual        = [];
        output.Slack       = [];
        output.problem     = -4;
        output.infostr     = yalmiperror(-4);
        output.solverinput = [];
        output.solveroutput= [];
        output.solvertime  = 0;
        return
    else
        spread = ub(integer_variables) - lb(integer_variables);
        bits = ceil(log2(spread+1));
        H = full(sparse(find(ismember(1:length(c),binary_variables)),1:length(binary_variables),ones(1,length(binary_variables)),length(c),length(binary_variables)));
        for i = 1:length(integer_variables)
            j = integer_variables(i);
            H(j,end+1:end+bits(i)) = 2.^(0:bits(i)-1);
        end
    end
    xbase = zeros(length(c),1);
    xbase(integer_variables) = lb(integer_variables);
    Aold = A;
    A = Aold*H;
    b = b-Aold*xbase;
    beq = beq-Aeq*xbase;
    Aeq = Aeq*H;
    c = H'*c;
else
    H = 1;
    xbase  = 0;
end

switch options.verbose
case 0
    options.bintprog.Display = 'off';
case 1
    options.bintprog.Display = 'final';
otherwise
    options.bintprog.Display = 'iter';
end

c = full(c);
A = full(A);
b = full(b);
Aeq = full(Aeq);
beq = full(beq);
ops = options.bintprog;
if options.savedebug
    save bintprogdebug c A b Aeq beq ops
end

[x,fmin,flag,output] = bintprog(c, A, b, Aeq, beq, x0,ops);

% Go back to integer variables
if ~isempty(x)
    x = xbase + H*x;
end

solvertime = etime(clock,solvertime);
problem = 0;

% No duals
D_struc = [];

% Check, currently not exhaustive...
switch flag
case 1
    problem = 0;
case {0,-4,-5,-6}
    problem = 3;
case -2
    problem = 1;
otherwise
    problem = -1;
end

infostr = yalmiperror(problem,'BINTPROG');       

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.Aeq = Aq;
    solverinput.beq = beq;
    solverinput.c = c;
    solverinput.options = options.linprog;
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
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;