function output = callgpposy(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
extended_variables = interfacedata.extended_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
mt      = interfacedata.monomtable;
variabletype = interfacedata.variabletype;

% *********************************
% What type of variables do we have
% *********************************
if isempty(variabletype)
    linear_variables = find((sum(abs(mt),2)==1) & (any(mt==1,2)));
else
    linear_variables = find(variabletype == 0);
end

% Same for fmincon,mosek,gpposy (for gpposy, we do not add bound constraints)
[prob,problem] = yalmip2geometric(options,F_struc,c,Q,K,[],[],mt,linear_variables,extended_variables);

if problem == 0
%     
    fixed = find(ub(linear_variables) == lb(linear_variables));
    if ~isempty(fixed)
        prob.G = [prob.G;-sparse(1:length(fixed),fixed,1,length(fixed),length(linear_variables))];
        prob.h = [prob.h;lb(linear_variables(fixed))];
    end
    ub(linear_variables(fixed)) = inf;
    lb(linear_variables(fixed)) = -inf;

    % Account for numerical problems in gpposy
    if ~isempty(lb);
        lb=lb(linear_variables);
        lb(lb<0)=1e-100;
        lb(lb==0) = 1e-100;
    end
    if ~isempty(ub);
        ub=ub(linear_variables);
        ub(isinf(ub))=1e100;
    end

    % Convert to gpposy
    A=prob.A;
    b=prob.b;
    G = prob.G;
    h = prob.h;
    szs=[];
    for i=0:max(prob.map)
        szs=[szs;nnz(find(i==prob.map))];
    end
    
    if szs(1) == 0
        % Feasibility problem not supported by GPPOSY
        % Just minimize sum of all variables
        A = [eye(size(A,2));A];
        b = [ones(size(A,2),1);b];
        szs(1) = size(A,2);
    end

    if options.savedebug
        save gpposydebug A b szs
    end
       
    solvertime = tic;
    [x,status,lambda,nu] = gpposy(A,b,szs,G,h,lb,ub,double(options.verbose)==0);
    solvertime = toc(solvertime);

    Primal = zeros(length(c),1);

    % Check, currently not exhaustive...
    switch lower(status)
        case 'solved'
            problem = 0;
            Primal(linear_variables) = x;
        case 'infeasible'
            problem = 1;
        case 'failed'
            problem = 4;
        otherwise
            problem = 9;
    end
else
    Primal = [];
    solvertime = [];
end


% Internal format for duals
Dual = [];

if options.savesolverinput
    solverinput.A   = A;
    solverinput.b   = b;
    solverinput.szs = szs;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.status = status;
else
    solveroutput = [];
end

% Standard interface
output = createOutputStructure(Primal,Dual,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);