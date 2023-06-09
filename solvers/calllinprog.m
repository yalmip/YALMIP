function output = calllinprog(interfacedata)

% Standard input interface
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;
x0      = interfacedata.x0;

showprogress('Calling LINPROG',options.showprogress);

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

switch options.verbose
case 0
    options.linprog.Display = 'off';
case 1
    options.linprog.Display = 'final';
otherwise
    options.linprog.Display = 'iter';
end

if isfield(options.linprog,'LargeScale')
    if ~isequal(options.linprog.LargeScale,'on')       
        c = full(c);
        A = full(A);
        b = full(b);
        Aeq = full(Aeq);
        beq = full(beq);
    end
end

if ~options.warmstart
    x0 = [];
end

if options.savedebug
    ops = options.linprog;
    save linprogdebug c A b Aeq beq lb ub ops x0
end

solvertime = tic;
[x,fmin,flag,output,lambda] = linprog(c, A, b, Aeq, beq, lb, ub, x0,options.linprog);
solvertime = toc(solvertime);
problem = 0;

% Internal format for duals
if isempty(lambda)
    D_struc = [];
else
    D_struc = [lambda.eqlin;lambda.ineqlin];
end
if isempty(x)
    x = zeros(length(c),1);
end

% Check, currently not exhaustive...
if flag==0
    problem = 3;
elseif flag == -2 | flag==-5
    problem = 1;
elseif flag == -3
    problem = 2;
else
    if flag>0
        problem = 0;
    else 
        if any((A*x-b)>sqrt(eps)) | any( abs(Aeq*x-beq)>sqrt(eps))
            problem = 1; % Likely to be infeasible
        else
            if c'*x<-1e10 % Likely unbounded
                problem = 2;
            else          % Probably convergence issues
                problem = 5;
            end
        end
    end
end

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.Aeq = Aeq;
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
output = createOutputStructure(x(:),D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);
