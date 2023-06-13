function output = callcopt(interfacedata)

options = interfacedata.options;
model = yalmip2copt(interfacedata);

if interfacedata.options.savedebug
    save coptdebug problem
end

if options.showprogress; showprogress('Call COPT', options.showprogress); end
solvertime = tic;
solution   = copt_solve(model, model.params);
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

        dims = model.conedata.K;

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
        problem = 0;
    case 'infeasible'
        problem = 1;
    case 'unbounded'
        problem = 2;
    case {'timeout', 'nodelimit'}
        problem = 3;
    case 'numerical'
        problem = 4;
    case 'interrupted'
        problem = 16;
    case 'unstarted'
        problem = -4;
    otherwise
        problem = -1;
end

if interfacedata.options.savesolverinput
    solverinput.model = model;
else
    solverinput = [];
end

if interfacedata.options.savesolveroutput
    solveroutput.result = solution;
else
    solveroutput = [];
end

output = createOutputStructure(x,D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);
