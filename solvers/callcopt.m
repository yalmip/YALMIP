function output = callcopt(interfacedata)

options = interfacedata.options;
problem = yalmip2copt(interfacedata);

if interfacedata.options.savedebug
    save coptdebug problem
end

if options.showprogress; showprogress('Call COPT', options.showprogress); end
solvertime = tic;
solution   = copt_solve(problem, problem.params);
solvertime = toc(solvertime);

if isfield(solution, 'rowmap')
    nrow = length(solution.rowmap);
    x = zeros(nrow, 1);
    if isfield(solution, 'psdpi')
        for i = 1:nrow
            if solution.rowmap(i) < 0
                x(i) = -solution.psdpi(-solution.rowmap(i));
            else
                x(i) = -solution.pi(solution.rowmap(i));
            end
        end
    end

    if isfield(solution, 'psdx')
        if isfield(solution, 'x')
            y = [solution.x; solution.psdx];
        else
            y = solution.psdx;
        end

        dims = problem.conedata.K;

        D_struc = y(1:dims.f + dims.l + sum(dims.q) + sum(dims.r));
        if isfield(solution, 'psdx')
            top = 1 + dims.f + dims.l + sum(dims.q) + sum(dims.r);
            for i = 1:length(dims.s)
                n = dims.s(i);
                sdpdual = y(top:top + n * (n + 1) / 2 - 1);
                Z = zeros(n);
                ind = find(tril(ones(n)));
                Z(ind) = sdpdual;
                Z = Z + Z';
                ind = find(speye(n));
                Z(ind) = Z(ind) / 2.0;
                D_struc = [D_struc; Z(:)];
                top = top + n * (n + 1) / 2;
            end
        end
    else
        D_struc = [];
    end

    dualstatus = solution.status;
    switch dualstatus
        case 'infeasible'
            solution.status = 'unbounded';
        case 'unbounded'
            solution.status = 'infeasible';
        otherwise
            solution.status = dualstatus;
    end
else
    numvars = length(interfacedata.c);

    if isfield(solution, 'x')
        x = solution.x(1:numvars);
    else
        x = zeros(numvars, 1);
    end

    if isfield(solution, 'pi')
        D_struc = -solution.pi;
    else
        D_struc = [];
    end
end

switch solution.status
    case 'optimal'
        status = 0;
    case 'infeasible'
        status = 1;
    case 'unbounded'
        status = 2;
    case {'timeout', 'nodelimit'}
        status = 3;
    case 'numerical'
        status = 4;
    case 'interrupted'
        status = 16;
    case 'unstarted'
        status = -4;
    otherwise
        status = -1;
end

if interfacedata.options.savesolverinput
    solverinput.model = problem;
else
    solverinput = [];
end

if interfacedata.options.savesolveroutput
    solveroutput.result = solution;
else
    solveroutput = [];
end

infostr = yalmiperror(status, interfacedata.solver.tag);

output  = createOutputStructure(x, D_struc, [], status, infostr, solverinput, solveroutput, solvertime);
