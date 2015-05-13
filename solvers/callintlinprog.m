function output = callintlinprog(interfacedata)

% Standard input interface
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;
x0      = interfacedata.x0;

showprogress('Calling INTLINPROG',options.showprogress);

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

intcon = union(interfacedata.integer_variables,interfacedata.binary_variables);
if ~isempty(interfacedata.binary_variables)
    lb(interfacedata.binary_variables) = max(0,lb(interfacedata.binary_variables));
    ub(interfacedata.binary_variables) = min(1,ub(interfacedata.binary_variables));
end
ops = options.intlinprog;

if options.verbose == 0
    ops.Display = 'off';
end

if options.savedebug    
    save intlinprogdebug c intcon A b Aeq beq lb ub ops x0
end
 
solvertime = tic;
[x,fval,exitflag,output] = intlinprog(c, intcon, A, b, Aeq, beq, lb, ub,ops);
solvertime = toc(solvertime);
problem = 0;

if isempty(x)
    x = zeros(length(c),1);
end

% Check, currently not exhaustive...
switch exitflag
    case 1
        problem = 0;
    case -2
        problem = 1;
    case -3
        problem = 2;       
    otherwise
        problem = 9;        
end
infostr = yalmiperror(problem,'INTLINPROG');       

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.Aeq = Aeq;
    solverinput.beq = beq;
    solverinput.c = c;
    solverinput.intcon = intcon;
    solverinput.options = ops;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
    solveroutput.output=output;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x(:),[],[],problem,infostr,solverinput,solveroutput,solvertime);