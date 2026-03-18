function output = call_highsmex(interfacedata)

options = interfacedata.options;
model = yalmip2quadprog(interfacedata);

% callhighs expects dense double vectors/matrices (or empty arrays)
c = full(double(model.c(:)));

if options.savedebug
    save highsmexdebug model options
end

if options.showprogress
    showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);
end

n = length(model.c);

% Build HiGHS row-wise bound form: L <= A*x <= U
if isempty(model.Aeq) && isempty(model.A)
    A = [];
    L = [];
    U = [];
else
    A = full(double([model.Aeq;model.A]));
    L = full(double([model.beq(:); -inf(length(model.b),1)]));
    U = full(double([model.beq(:); model.b(:)]));
end

Q = [];
if nnz(model.Q) > 0
    Q = full(double(model.Q));
end

lb = full(double(model.lb(:)));
ub = full(double(model.ub(:)));

isMixedInteger = ~isempty(interfacedata.integer_variables) || ~isempty(interfacedata.binary_variables);
if isMixedInteger && nnz(model.Q) > 0
    problem = -4;
    infostr = yalmiperror(problem,interfacedata.solver.tag);
    if options.savesolverinput
        solverinput.model = model;
        solverinput.c = c;
        solverinput.A = A;
        solverinput.L = L;
        solverinput.U = U;
        solverinput.lb = lb;
        solverinput.ub = ub;
        solverinput.Q = Q;
    else
        solverinput = [];
    end
    if options.savesolveroutput
        solveroutput.info.model_status_string = 'MIQP not supported by HiGHS';
    else
        solveroutput = [];
    end
    output = createOutputStructure(nan(n,1),[],[],problem,infostr,solverinput,solveroutput,0);
    return
end

integrality = repmat({'c'},n,1);
if ~isempty(interfacedata.integer_variables)
    integrality(interfacedata.integer_variables) = {'i'};
end
if ~isempty(interfacedata.binary_variables)
    integrality(interfacedata.binary_variables) = {'i'};
    lb(interfacedata.binary_variables) = max(0,lb(interfacedata.binary_variables));
    ub(interfacedata.binary_variables) = min(1,ub(interfacedata.binary_variables));
end

allcontinuous = all(strcmp(integrality,'c'));
if allcontinuous
    integralityArg = [];
else
    if exist('string','builtin') || exist('string','file')
        integralityArg = string(integrality);
    else
        integralityArg = integrality;
    end
end

highsops = options.highsmex;
if isempty(highsops)
    highsops = struct;
end
if isstruct(highsops) && options.verbose == 0
    highsops.log_to_console = false;
elseif isstruct(highsops) && options.verbose > 0
    highsops.log_to_console = true;
end

setSoln = [];
if options.usex0 && ~isempty(model.x0) && all(isfinite(model.x0))
    setSoln = model.x0(:);
end

solvertime = tic;
problem = 0;
soln = struct('col_value',nan(n,1),'row_dual',[]);
info = struct('model_status_string','Solve error','run_time',0);
optsused = [];
basis = [];

[soln,info,optsused,basis] = callhighs(c,A,L,U,lb,ub,Q,integralityArg,highsops,[],setSoln,[]);
if isfield(info,'run_time') && ~isempty(info.run_time)
    solvertime = info.run_time;
else
    solvertime = toc(solvertime);
end
problem = map_highsmex_status(info);

if isfield(soln,'col_value') && ~isempty(soln.col_value)
    Primal = soln.col_value(:);
else
    Primal = nan(n,1);
end

Dual = [];
if options.saveduals && allcontinuous
    if isfield(soln,'row_dual') && ~isempty(soln.row_dual)
        % HiGHS row dual sign is opposite to YALMIP's Ax <= b convention.
        Dual = -soln.row_dual(:);
    end
end

infostr = yalmiperror(problem,interfacedata.solver.tag);

if options.savesolverinput
    solverinput.model = model;
    solverinput.c = c;
    solverinput.A = A;
    solverinput.L = L;
    solverinput.U = U;
    solverinput.lb = lb;
    solverinput.ub = ub;
    solverinput.Q = Q;
    solverinput.integrality = integralityArg;
    solverinput.options = highsops;
    solverinput.setSoln = setSoln;
else
    solverinput = [];
end

if options.savesolveroutput
    solveroutput.soln = soln;
    solveroutput.info = info;
    solveroutput.options = optsused;
    solveroutput.basis = basis;
else
    solveroutput = [];
end

output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);


function problem = map_highsmex_status(info)

problem = -1;
if ~isstruct(info) || ~isfield(info,'model_status_string') || isempty(info.model_status_string)
    return
end

status = lower(char(info.model_status_string));

if ~isempty(strfind(status,'optimal')) || ~isempty(strfind(status,'objective target'))
    problem = 0;
elseif ~isempty(strfind(status,'infeasible')) && ~isempty(strfind(status,'unbounded'))
    problem = 12;
elseif ~isempty(strfind(status,'infeasible'))
    problem = 1;
elseif ~isempty(strfind(status,'unbounded'))
    problem = 2;
elseif ~isempty(strfind(status,'time limit')) || ~isempty(strfind(status,'iteration limit')) || ~isempty(strfind(status,'solution limit')) || ~isempty(strfind(status,'objective bound'))
    problem = 3;
elseif ~isempty(strfind(status,'interrupt'))
    problem = 16;
elseif ~isempty(strfind(status,'error')) || ~isempty(strfind(status,'notset'))
    problem = 9;
end
