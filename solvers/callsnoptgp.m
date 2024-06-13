function output = callsnoptgp(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
extended_variables = interfacedata.extended_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
mt      = interfacedata.monomtable;

% *********************************
% What type of variables do we have
% *********************************
linear_variables = find((sum(abs(mt),2)==1) & (any(mt==1,2)));
nonlinear_variables = setdiff((1:size(mt,1))',linear_variables);
sigmonial_variables = find(any(0>mt,2) | any(mt-fix(mt),2));

% Convert to common format for fmincon, mosek and gpposy
ubtemp = ub;
lbtemp = lb;
fixed = find(ub(linear_variables) == lb(linear_variables));
ubtemp(linear_variables(fixed)) = inf;
lbtemp(linear_variables(fixed)) = -inf;
[prob,problem] = yalmip2geometric(options,F_struc,c,Q,K,ubtemp,lbtemp,mt,linear_variables,extended_variables);

%Add equalities for fixed variables
fixed_in_bound = [];
fixed_bound=[];
if ~isempty(fixed)
    if 0
        % Remove variables
        fixfactors = prod(repmat(ub(linear_variables(fixed(:)))',size(prob.A,1),1).^prob.A(:,fixed),2);
        prob.A(:,fixed) = [];
        prob.b = (prob.b) .* fixfactors;
        for j = 1:max(prob.map)
            in_this = find(prob.map == j);
            if ~isempty(in_this)
                constants = find(~any(prob.A(in_this,:),2));
                if ~isempty(constants)
                    prob.b(in_this) = prob.b(in_this)./(1 - sum(prob.b(in_this(constants))));
                end
            end
        end
        fixed_in_bound = linear_variables(fixed);
        fixed_bound = ub(linear_variables(fixed));
        linear_variables(fixed) = [];
        em = find(~any(prob.A,2));
        prob.A(em,:) = [];
        prob.b(em,:) = [];
        nn = max(prob.map);
        prob.map(em,:) = [];
        length(unique(prob.map)) - (nn + 1)
        1;
    else
        prob.G = [prob.G;-sparse(1:length(fixed),fixed,1,length(fixed),length(linear_variables))];
        prob.h = [prob.h;lb(linear_variables(fixed))];
    end
end

%something failed, perhaps a QP
if problem & isempty(sigmonial_variables)
    % This is an LP or QP!
    % Go to standard fmincon
    if options.verbose
        disp('Conversion to geometric program failed. Trying general non-convex model in fmincon');
        disp(' ');
    end
    interfacedata.solver.tag = strrep(interfacedata.solver.tag,'-geometric','');
    output = callsnopt(interfacedata);
    return
end

if problem == 0
    if isempty(x0)
        x0 = zeros(length(linear_variables),1);
    else
        x0 = x0(linear_variables);
    end

    % Fake logarithm (extend linearly for small values)
    ind = find(x0<1e-2);
    x0(ind) = exp(log(1e-2)+(x0(ind)-1e-2)/1e-2);
    x0 = log(x0);

    % Clean up the bounds (from branch and bound)
    % Note, these bounds are in the
    % logarithmic variables.
    if ~isempty(lb)
        lb = lb(linear_variables);
        ind = find(lb<1e-2);
        lb(ind) = exp(log(1e-2)+(lb(ind)-1e-2)/1e-2);
        lb = log(lb+sqrt(eps));
    end
    if ~isempty(ub)
        ub = ub(linear_variables);
        ind = find(ub<1e-2);
        ub(ind) = exp(log(1e-2)+(ub(ind)-1e-2)/1e-2);
        ub = log(ub+sqrt(eps));
    end

    if options.savedebug
        ops = options.fmincon;
        save snoptgpdebug prob x0 ops lb ub
    end

    prob = precalcgpstruct(prob);

    solveopt = 1;
    xlow = lb;
    xupp = ub;
    xmul   = zeros(length(xupp),1);
    xstate = zeros(length(xupp),1);

    Fupp = [ inf;
        repmat(0,max(prob.map),1);
        repmat(0,size(prob.G,1),1)];

    Flow = [ -inf;
        repmat(-inf,max(prob.map),1);
        repmat(0,size(prob.G,1),1)];
    
    Fmul   = zeros(length(Fupp),1);
    Fstate = zeros(length(Fupp),1);
    ObjAdd = 0;
    ObjRow = 1;
    A = [];
    iAfun = [];
    jAvar = [];
    G = ones(size(Fupp,1),length(xupp));
    [iGfun,jGvar] = find(G);
    usrf = 'snoptgp_callback';
    snoptgp_callback([],prob);
    solvertime = tic;
    if options.verbose == 0
        evalc('[xout,F,xmul,Fmul,inform, xstate, Fstate, ns, ninf, sinf, mincw, miniw, minrw] = snoptcmex( solveopt, x0, xlow, xupp, xmul, xstate, Flow, Fupp, Fmul, Fstate,ObjAdd, ObjRow, A, iAfun(:), jAvar(:),iGfun(:), jGvar(:), usrf );');
    else
        [xout,F,xmul,Fmul,inform, xstate, Fstate, ns, ninf, sinf, mincw, miniw, minrw] = snoptcmex( solveopt, x0, xlow, xupp, xmul, xstate, Flow, Fupp, Fmul, Fstate,ObjAdd, ObjRow, A, iAfun(:), jAvar(:),iGfun(:), jGvar(:), usrf );
    end
    solvertime = toc(solvertime);
    lambda = Fmul(2:end);    

    x = zeros(length(c),1);
    x(linear_variables) = exp(xout);
    x(fixed_in_bound) = fixed_bound;

    problem = 0;

    % Internal format for duals
    D_struc = [];

    % Check, currently not exhaustive...
    switch inform
        case {0,1,2,3}
            problem = 0;
        case {11,12,13,14}
            problem = 1;
        otherwise
    end

else
    x = [];
    solvertime = 0;
end

% Internal format for duals (currently not supported in GP)
D_struc = [];

if options.savesolverinput
    solverinput.A = [];
    solverinput.b = [];
    solverinput.Aeq = [];
    solverinput.beq = [];
    solverinput.options = options.fmincon;
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
output = createOutputStructure(x,D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);