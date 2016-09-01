function output = bnb(p)
%BNB          General branch-and-bound scheme for conic programs
%
% BNB applies a branch-and-bound scheme to solve mixed integer
% conic programs (LP, QP, SOCP, SDP) and mixed integer geometric programs.
%
% BNB is never called by the user directly, but is called by
% YALMIP from SOLVESDP, by choosing the solver tag 'bnb' in sdpsettings.
%
% BNB is used if no other mixed integer solver is found, and
% is only useful for very small problems, due to its simple
% and naive implementation.
%
% The behaviour of BNB can be altered using the fields
% in the field 'bnb' in SDPSETTINGS
%
% bnb.branchrule   Deceides on what variable to branch
%                   'max'     : Variable furthest away from being integer
%                   'min'     : Variable closest to be being integer
%                   'first'   : First variable (lowest variable index in YALMIP)
%                   'last'    : Last variable (highest variable index in YALMIP)
%                   'weight'  : See manual
%
% bnb.method       Branching strategy
%                   'depth'   : Depth first
%                   'breadth' : Breadth first
%                   'best'    : Expand branch with lowest lower bound
%                   'depthX'  : Depth until integer solution found, then X (e.g 'depthbest')
%
% solver           Solver for the relaxed problems (standard solver tag, see SDPSETTINGS)
%
% maxiter          Maximum number of nodes explored
%
% inttol           Tolerance for declaring a variable as integer
%
% feastol          Tolerance for declaring constraints as feasible
%
% gaptol           Exit when (upper bound-lower bound)/(1e-3+abs(lower bound)) < gaptol
%
% round            Round variables smaller than bnb.inttol
%
%
% See also SOLVESDP, BINVAR, INTVAR, BINARY, INTEGER

% ********************************
%% INITIALIZE DIAGNOSTICS IN YALMIP
% ********************************
bnbsolvertime = clock;
showprogress('Branch and bound started',p.options.showprogress);

% ********************************
%% We might have a GP : pre-calc
% ********************************
p.nonlinear = find(~(sum(p.monomtable~=0,2)==1 & sum(p.monomtable,2)==1));
p.nonlinear = union(p.nonlinear,p.evalVariables);

% ********************************
% This field is only used in bmibnb, which uses the same sub-functions as
% bnb
% ********************************
p.high_monom_model = [];

% ********************************
%% Define infinite bounds
% ********************************
if isempty(p.ub)
    p.ub = repmat(inf,length(p.c),1);
end
if isempty(p.lb)
    p.lb = repmat(-inf,length(p.c),1);
end

% ********************************
%% Extract bounds from model
% ********************************
if ~isempty(p.F_struc)
    [lb,ub,used_rows_eq,used_rows_lp] = findulb(p.F_struc,p.K);
    if ~isempty(used_rows_lp)
        used_rows_lp = used_rows_lp(~any(full(p.F_struc(p.K.f + used_rows_lp,1+p.nonlinear)),2));
        if ~isempty(used_rows_lp)
            lower_defined = find(~isinf(lb));
            if ~isempty(lower_defined)
                p.lb(lower_defined) = max(p.lb(lower_defined),lb(lower_defined));
            end
            upper_defined = find(~isinf(ub));
            if ~isempty(upper_defined)
                p.ub(upper_defined) = min(p.ub(upper_defined),ub(upper_defined));
            end
            p.F_struc(p.K.f + used_rows_lp,:)=[];
            p.K.l = p.K.l - length(used_rows_lp);
        end
    end
    
    if ~isempty(used_rows_eq)
        used_rows_eq = used_rows_eq(~any(full(p.F_struc(used_rows_eq,1+p.nonlinear)),2));
        if ~isempty(used_rows_eq)
            lower_defined = find(~isinf(lb));
            if ~isempty(lower_defined)
                p.lb(lower_defined) = max(p.lb(lower_defined),lb(lower_defined));
            end
            upper_defined = find(~isinf(ub));
            if ~isempty(upper_defined)
                p.ub(upper_defined) = min(p.ub(upper_defined),ub(upper_defined));
            end
            p.F_struc(used_rows_eq,:)=[];
            p.K.f = p.K.f - length(used_rows_eq);
        end
    end
end

% ********************************
%% ADD CONSTRAINTS 0<x<1 FOR BINARY
% ********************************
if ~isempty(p.binary_variables)
    p.ub(p.binary_variables) =  min(p.ub(p.binary_variables),1);
    p.lb(p.binary_variables) =  max(p.lb(p.binary_variables),0);
    
   % godown = find(p.ub(p.binary_variables) < 1);
   % goup   = find(p.lb(p.binary_variables) > 0);
   % p.ub(p.binary_variables(godown)) = 0;
   % p.lb(p.binary_variables(goup)) = 1;
end

%p.lb(p.integer_variables) = ceil(p.lb(p.integer_variables));
%p.ub(p.integer_variables) = floor(p.ub(p.integer_variables));
p = update_integer_bounds(p);

if ~isempty(p.semicont_variables)
    redundant = find(p.lb<=0 & p.ub>=0);
    p.semicont_variables = setdiff(p.semicont_variables,redundant);
    % Now relax the model and generate hull including 0
    p.semibounds.lb = p.lb(p.semicont_variables);
    p.semibounds.ub = p.ub(p.semicont_variables);
    p.lb(p.semicont_variables) = min(p.lb(p.semicont_variables),0);
    p.ub(p.semicont_variables) = max(p.ub(p.semicont_variables),0);
end

% Could be some nonlinear terms (although these problems are recommended to
% be solved using BMIBNB
p = compile_nonlinear_table(p);
p = updatemonomialbounds(p);

% *******************************
%% PRE-SOLVE (nothing fancy coded)
% *******************************
pss=[];
p = propagate_bounds_from_equalities(p);

if p.K.f > 0
    pp = p;
    r = find(p.lb == p.ub);
    pp.F_struc(:,1) = pp.F_struc(:,1) + pp.F_struc(:,r+1)*p.lb(r);
    pp.F_struc(:,r+1)=[];
    pp.lb(r)=[];
    pp.ub(r)=[];
    pp.variabletype(r)=[];
    % FIXME: This is lazy, should update new list
    pp.binary_variables = [];
    pp.integer_variables = [];
    pp = propagate_bounds_from_equalities(pp);
    other = setdiff(1:length(p.lb),r);
    p.lb(other) = pp.lb;
    p.ub(other) = pp.ub;
    p = update_integer_bounds(p);
    redundant = find(~any(pp.F_struc(1:p.K.f,2:end),2));
    if any(p.F_struc(redundant,1)<0)
        p.feasible = 0;
    else
        p.F_struc(redundant,:)=[];
        p.K.f = p.K.f - length(redundant);
    end
end

if isempty(p.nonlinear)
    if p.K.f>0
        Aeq = -p.F_struc(1:p.K.f,2:end);
        beq = p.F_struc(1:p.K.f,1);
        A = [Aeq;-Aeq];
        b = [beq;-beq];
        [p.lb,p.ub,redundant,pss] = tightenbounds(A,b,p.lb,p.ub,p.integer_variables,p.binary_variables,ones(length(p.lb),1));
    end
    pss=[];
    if p.K.l>0
        A = -p.F_struc(1+p.K.f:p.K.f+p.K.l,2:end);
        b = p.F_struc(1+p.K.f:p.K.f+p.K.l,1);
        [p.lb,p.ub,redundant,pss] = tightenbounds(A,b,p.lb,p.ub,p.integer_variables,p.binary_variables,ones(length(p.lb),1));
        if length(redundant)>0
            pss.AL0A(redundant,:)=[];
            pss.AG0A(redundant,:)=[];
            p.F_struc(p.K.f+redundant,:)=[];
            p.K.l = p.K.l - length(redundant);
        end
    end
end

% Silly redundancy
p = updatemonomialbounds(p);
p = propagate_bounds_from_equalities(p);
if p.K.l > 0
    b = p.F_struc(1+p.K.f:p.K.l+p.K.f,1);
    A = -p.F_struc(1+p.K.f:p.K.l+p.K.f,2:end);
    redundant = find(((A>0).*A*(p.ub-p.lb) - (b-A*p.lb) <= 0));
    if ~isempty(redundant)
        p.F_struc(p.K.f + redundant,:) = [];
        p.K.l = p.K.l - length(redundant);
    end
end

% *******************************
%% PERTURBATION OF LINEAR COST
% *******************************
p.corig = p.c;
if nnz(p.Q)==0 & isequal(p.K.m,0)
%     g = randn('seed');
%     randn('state',1253); %For my testing, I keep this the same...
%     % This perturbation has to be better. Crucial for many real LP problems
%     p.c = (p.c).*(1+randn(length(p.c),1)*1e-4);
%     randn('seed',g);
end

% *******************************
%% Display logics
% 0 : Silent
% 1 : Display branching
% 2 : Display node solver prints
% *******************************
switch max(min(p.options.verbose,3),0)
    case 0
        p.options.bnb.verbose = 0;
    case 1
        p.options.bnb.verbose = 1;
        p.options.verbose = 0;
    case 2
        p.options.bnb.verbose = 2;
        p.options.verbose = 0;
    case 3
        p.options.bnb.verbose = 2;
        p.options.verbose = 1;
    otherwise
        p.options.bnb.verbose = 0;
        p.options.verbose = 0;
end

% *******************************
%% Figure out the weights if any
% *******************************
try % Probably buggy first version...
    if ~isempty(p.options.bnb.weight)
        weightvar = p.options.bnb.weight;
        if isa(weightvar,'sdpvar')
            if (prod(size(weightvar)) == 1)
                weight = ones(length(p.c),1);
                for i = 1:length(p.c)
                    weight(i,1) = full(getbasematrix(weightvar,p.used_variables(i)));
                end
                p.weight = weight;
            else
                error('Weight should be an SDPVAR scalar');
            end
        else
            error('Weight should be an SDPVAR scalar');
        end
    else
        p.weight = ones(length(p.c),1);
    end
catch
    disp('Something wrong with weights. Please report bug');
    p.weight = ones(length(p.c),1);
end

% *******************************
%% START BRANCHING
% *******************************
setuptime = etime(clock,bnbsolvertime);
bnbsolvertime = clock;
[x_min,solved_nodes,lower,upper,profile,diagnostics] = branch_and_bound(p,pss);
bnbsolvertime =  etime(clock,bnbsolvertime);
output.solvertime   = setuptime + bnbsolvertime;

% **********************************
%% CREATE SOLUTION
% **********************************
if diagnostics == -4
    output.problem = -4;
else
    output.problem = 0;
    if isinf(upper)
        output.problem = 1;
    end
    if isinf(-lower)
        output.problem = 2;
    end
    if solved_nodes == p.options.bnb.maxiter
        output.problem = 3;
    end
end
output.solved_nodes = solved_nodes;
output.Primal      = x_min;
output.Dual        = [];
output.Slack       = [];
if  output.problem == -4
    output.infostr      = yalmiperror(output.problem,[p.solver.lower.tag '-' p.solver.lower.version]);
else
    output.infostr      = yalmiperror(output.problem,'BNB');
end
output.solverinput  = 0;
if p.options.savesolveroutput
    output.solveroutput.setuptime = setuptime;
    output.solveroutput.localsolvertime = profile.local_solver_time;
    output.solveroutput.branchingtime = bnbsolvertime;
    output.solveroutput.solved_nodes = solved_nodes;
    output.solveroutput.lower = lower;
    output.solveroutput.upper = upper;
else
    output.solveroutput =[];
end
%% --

function [x_min,solved_nodes,lower,upper,profile,diagnostics] = branch_and_bound(p,pss)

% *******************************
% We don't need this
% *******************************
p.options.savesolveroutput = 0;
p.options.saveduals = 0;
p.options.dimacs = 0;
diagnostics = 0;
% *******************************
% Tracking performance etc
% *******************************
profile.local_solver_time = 0;

% *************************************************************************
% We save this to re-use some stuff in fmincon
% *************************************************************************
p.options.savesolverinput = 1;

% *******************************
%% SET-UP ROOT PROBLEM
% *******************************
p.depth = 0;
p.lower = NaN;
% Does the user want to create his own initial guess
if p.options.usex0
    [x_min,upper] = initializesolution(p);
    if isinf(upper)
        % Try to initialize to lowerbound+ upperbound. fmincon really
        % doesn't like zero initial guess, despite having bounds available
        x_min   = zeros(length(p.c),1);
        violates_finite_bounds = ((x_min < p.lb) | (x_min < p.ub));
        violates_finite_bounds = find(violates_finite_bounds & ~isinf(p.lb) & ~isinf(p.ub));
        x_min(violates_finite_bounds) = (p.lb(violates_finite_bounds) + p.ub(violates_finite_bounds))/2;
        x_min  = setnonlinearvariables(p,x_min);
    end
    p.x0    = x_min;
else
    upper   = inf;
    x_min   = zeros(length(p.c),1);
    violates_finite_bounds = ((x_min < p.lb) | (x_min < p.ub));
    violates_finite_bounds = find(violates_finite_bounds & ~isinf(p.lb) & ~isinf(p.ub));
    x_min(violates_finite_bounds) = (p.lb(violates_finite_bounds) + p.ub(violates_finite_bounds))/2;
    x_min  = setnonlinearvariables(p,x_min);
    p.x0    = x_min;
end


% *******************************
%% Global stuff
% *******************************
lower   = NaN;
stack   = [];

% *******************************
%% Create function handle to solver
% *******************************
lowersolver = p.solver.lower.call;
uppersolver = p.options.bnb.uppersolver;

% *******************************
%% INVARIANT PROBLEM DATA
% *******************************
c = p.corig;
Q = p.Q;
f = p.f;
integer_variables = p.integer_variables;
solved_nodes = 0;

semicont_variables = p.semicont_variables;

gap = inf;
node = 1;

if p.options.bnb.presolve
    savec = p.c;
    saveQ = p.Q;
    p.Q = p.Q*0;
    
    n = length(p.c);
    saveBinary = p.binary_variables;
    saveInteger = p.integer_variables;
    p.binary_variables = [];
    p.integer_variables = [];;
    
    for i = 1:length(c)
        p.c = eyev(n,i);
        output = feval(lowersolver,p);
        if output.problem == 0
            p.lb(i) = max(p.lb(i),output.Primal(i));
        end
        p.c = -eyev(n,i);
        output = feval(lowersolver,p);
        if output.problem == 0
            p.ub(i) = min(p.ub(i),output.Primal(i));
        end
        p.lb(saveBinary) = ceil(p.lb(saveBinary)-1e-3);
        p.ub(saveBinary) = floor(p.ub(saveBinary)+1e-3);
    end
    p.binary_variables = saveBinary;
    p.integer_variables = saveInteger;
    
    p.Q = saveQ;
    p.c = savec;
end

% ************************************************
% Some hacks to speed up solver calls
% Only track solver-time if user wants profile
% ************************************************
p.getsolvertime = p.options.bnb.profile;

% *******************************
%% DISPLAY HEADER
% *******************************
originalDiscrete = [p.integer_variables(:);p.binary_variables(:)];
originalBinary   = p.binary_variables(:);

if nnz(Q)==0 & (nnz(p.c-fix(p.c))==0) & isequal(p.K.m,0)
    can_use_ceil_lower = all(ismember(find(p.c),originalDiscrete));
else
    can_use_ceil_lower = 0;
end

if p.options.bnb.verbose
    
    pc = p.problemclass;
    non_convex_obj = pc.objective.quadratic.nonconvex | pc.objective.polynomial;
    non_convex_constraint =  pc.constraint.equalities.quadratic | pc.constraint.inequalities.elementwise.quadratic.nonconvex;
    non_convex_constraint =  non_convex_constraint | pc.constraint.equalities.polynomial | pc.constraint.inequalities.elementwise.polynomial;
    
    possiblynonconvex = non_convex_obj | non_convex_constraint;
    if ~isequal(p.solver.lower.version,'')
        p.solver.lower.tag = [p.solver.lower.tag '-' p.solver.lower.version];
    end
    
    disp('* Starting YALMIP integer branch & bound.');
    disp(['* Lower solver   : ' p.solver.lower.tag]);
    disp(['* Upper solver   : ' p.options.bnb.uppersolver]);
    disp(['* Max iterations : ' num2str(p.options.bnb.maxiter)]);
    
    if possiblynonconvex & p.options.warning
        disp(' ');
        disp('Warning : The continuous relaxation may be nonconvex. This means ');
        disp('that the branching process is not guaranteed to find a');
        disp('globally optimal solution, since the lower bound can be');
        disp('invalid. Hence, do not trust the bound or the gap...')
        
    end
end
if p.options.bnb.verbose;            disp(' Node       Upper       Gap(%)      Lower    Open');end;

if nnz(Q)==0 & nnz(c)==1 & isequal(p.K.m,0)
    p.simplecost = 1;
else
    p.simplecost = 0;
end

poriginal = p;
p.cuts = [];

%% MAIN LOOP
% p.options.rounding = [1 1 1];

if p.options.bnb.nodefix & (p.K.s(1)>0)
    top=1+p.K.f+p.K.l+sum(p.K.q);
    for i=1:length(p.K.s)
        n=p.K.s(i);
        for j=1:size(p.F_struc,2)-1;
            X=full(reshape(p.F_struc(top:top+n^2-1,j+1),p.K.s(i),p.K.s(i)));
            X=(X+X')/2;
            v=real(eig(X+sqrt(eps)*eye(length(X))));
            if all(v>=0)
                sdpmonotinicity(i,j)=-1;
            elseif all(v<=0)
                sdpmonotinicity(i,j)=1;
            else
                sdpmonotinicity(i,j)=nan;
            end
        end
        top=top+n^2;
    end
else
    sdpmonotinicity=[];
end

% Try to find sum(d_i) = 1
sosgroups = {};
sosvariables = [];
if p.K.f > 0 & ~isempty(p.binary_variables)
    nbin = length(p.binary_variables);
    Aeq = -p.F_struc(1:p.K.f,2:end);
    beq = p.F_struc(1:p.K.f,1);
    notbinary_var_index = setdiff(1:length(p.lb),p.binary_variables);
    only_binary = ~any(Aeq(:,notbinary_var_index),2);
    Aeq_bin = Aeq(find(only_binary),p.binary_variables);
    beq_bin = beq(find(only_binary),:);
    % Detect groups with constraints sum(d_i) == 1
    sosgroups = {};
    for i = 1:size(Aeq_bin,1)
        if beq_bin(i) == 1
            [ix,jx,sx] = find(Aeq_bin(i,:));
            if all(sx == 1)
                sosgroups{end+1} = p.binary_variables(jx);
                sosvariables = [sosvariables p.binary_variables(jx)];
            end
        end
    end
end

pid = 0;
lowerhist = [];
upperhist = [];
p.fixedvariable = [];
p.fixdir = '';
p.sosgroups = sosgroups;
p.sosvariables = sosvariables;
oldp = p;
while ~isempty(node) & (solved_nodes < p.options.bnb.maxiter) & (isinf(lower) | gap>p.options.bnb.gaptol)
    
    % ********************************************
    % Adjust variable bound based on upper bound
    % ********************************************
    % This code typically never runs but can be turned on
    % using options.bnb.nodetight and bnb.nodefix.
    if ~isinf(upper) & ~isnan(lower)
        [p,poriginal,stack] = pruneglobally(p,poriginal,upper,lower,stack,x);
        [p,poriginal,stack] = fixvariables(p,poriginal,upper,lower,stack,x_min,sdpmonotinicity);
        stack = prunecardinality(p,poriginal,stack,lower,upper);
    end
   
    % ********************************************
    % BINARY VARIABLES ARE FIXED ALONG THE PROCESS
    % ********************************************
    binary_variables  = p.binary_variables;
    
    % ********************************************
    % SO ARE SEMI VARIABLES
    % ********************************************
    semicont_variables  = p.semicont_variables;
    
    % ********************************************
    % ASSUME THAT WE WON'T FATHOME
    % ********************************************
    keep_digging = 1;
    message = '';
    
    
    % *************************************
    % SOLVE NODE PROBLEM
    % *************************************
    if any(p.ub<p.lb - 1e-12)
        x = zeros(length(p.c),1);
        output.Primal = x;
        output.problem=1;
    else
        p.x_min = x_min;
        relaxed_p = p;
        relaxed_p.integer_variables = [];
        relaxed_p.binary_variables = [];
        relaxed_p.semicont_variables = [];
        relaxed_p.ub(p.ub<p.lb) = relaxed_p.lb(p.ub<p.lb);
        if upper<inf & length(poriginal.binary_variables)==length(poriginal.c) & p.K.f == 0
            % Cut away current best
            % FIXME: Generalize
            positive = find(x_min==1);
            zero = find(x_min==0);
            % Add cut c'*x < c*xmin,
            cc = poriginal.c;
            cc(positive) = ceil(cc(positive));
            cc(zero) = floor(cc(zero));
            relaxed_p.K.l =  relaxed_p.K.l+1;
            relaxed_p.F_struc = [-1+sum(cc(positive)) -cc';relaxed_p.F_struc];            
        end
       
        output = bnb_solvelower(lowersolver,relaxed_p,upper+abs(upper)*1e-2+1e-4,lower);
                             
        if p.options.bnb.profile
            profile.local_solver_time  = profile.local_solver_time + output.solvertime;
        end
        
        % A bit crappy code to exploit computations that were done in the
        % call to fmincon...
        if isfield(output,'solverinput')
            if isfield(output.solverinput,'model')
                if isfield(output.solverinput.model,'fastdiff')
                    p.fastdiff = output.solverinput.model.fastdiff;
                end
            end
        end
        
        
        if output.problem == -4
            diagnostics = -4;
            x = nan+zeros(length(p.lb),1);
        else
            if isempty(output.Primal)
                output.Primal = zeros(length(p.c),1);
            end
            try
                x  = setnonlinearvariables(p,output.Primal);
            catch
                1
            end
 
            if(p.K.l>0) & any(p.F_struc(p.K.f+1:p.K.f+p.K.l,:)*[1;x]<-1e-5)
                output.problem = 1;
            elseif output.problem == 5 & ~checkfeasiblefast(p,x,p.options.bnb.feastol)
                output.problem = 1;
            end            
        end
    end
    
    solved_nodes = solved_nodes+1;
    
    % **************************************
    % THIS WILL BE INTIAL GUESS FOR CHILDREN
    % **************************************
    p.x0 = x;
    
    % *************************************
    % ANY INTEGERS? ROUND?
    % *************************************
    non_integer_binary = abs(x(binary_variables)-round(x(binary_variables)))>p.options.bnb.inttol;
    non_integer_integer = abs(x(integer_variables)-round(x(integer_variables)))>p.options.bnb.inttol;
    if p.options.bnb.round
        x(binary_variables(~non_integer_binary))   = round(x(binary_variables(~non_integer_binary)));
        x(integer_variables(~non_integer_integer)) = round(x(integer_variables(~non_integer_integer)));
    end
    non_integer_binary = find(non_integer_binary);
    non_integer_integer = find(non_integer_integer);
    if isempty(p.semicont_variables)
        non_semivar_semivar=[];
    else
        non_semivar_semivar = find(~(abs(x(p.semicont_variables))<p.options.bnb.inttol | (x(p.semicont_variables)>p.semibounds.lb & x(p.semicont_variables)<=p.semibounds.ub)));
    end
    
    try
        x  = setnonlinearvariables(p,x);
    catch
    end
    
    TotalIntegerInfeas = sum(abs(round(x(non_integer_integer))-x(non_integer_integer)));
    TotalBinaryInfeas = sum(abs(round(x(non_integer_binary))-x(non_integer_binary)));
    
    % *************************************
    % NODE HEURISTICS (NOTHING CODED)
    % *************************************
    should_be_tight = find([p.lb == p.ub]);
    if ~isempty(should_be_tight)
        % FIX for problems that only report numerical problems but violate
        % binary
        if max(abs(p.lb(should_be_tight)-x(should_be_tight)))>p.options.bnb.inttol
            output.problem = 1;
        end
    end
    if output.problem==0 | output.problem==3 | output.problem==4
        cost = computecost(f,c,Q,x,p);
        
        if output.problem~=1
            if isnan(lower)
                lower = cost;
            end
            
            if isfield(p.options,'plottruss')
                if p.options.plottruss
                    plottruss(1,'Relaxed node',poriginal,x);
                end
            end
            
            if cost <= upper & ~(isempty(non_integer_binary) & isempty(non_integer_integer) & isempty(non_semivar_semivar))
                poriginal.upper = upper;
                poriginal.lower = lower;
                [upper1,x_min1] = feval(uppersolver,poriginal,output,p);
                if upper1 < upper
                    
                    if isfield(p.options,'plottruss')
                        if p.options.plottruss
                            plottruss(3,'Best binary solution',poriginal,x_min1);
                        end
                    end
                    
                    x_min = x_min1;
                    upper = upper1;
                    [stack,stacklower] = prune(stack,upper,p.options,solved_nodes,p);
                    lower = min(lower,stacklower);
                    [p,poriginal,stack] = pruneglobally(p,poriginal,upper,lower,stack,x_min);
                    [p,poriginal,stack] = fixvariables(p,poriginal,upper,lower,stack,x_min,sdpmonotinicity);
                    
                end
            elseif isempty(non_integer_binary) && isempty(non_integer_integer) && isempty(non_semivar_semivar)
            end
        end
    end
    p = adaptivestrategy(p,upper,solved_nodes);
    
    % *************************************
    % CHECK FATHOMING POSSIBILITIES
    % *************************************
    feasible = 1;
    
    switch output.problem
        case -1;
            % Solver behaved weird. Make sure we continue digging
            keep_digging = 1;
            feasible = 1;
            cost = lower;
            x = p.lb + (p.ub-p.lb)*(1/pi);
        case 0
            if can_use_ceil_lower
                lower = ceil(lower-1e-8);
            end
        case {1,12,-4}
            keep_digging = 0;
            cost = inf;
            feasible = 0;
        case 2
            cost = -inf;
        otherwise
            % This part has to be much more robust
            cost = f+c'*x+x'*Q*x;
    end
    
    
    
    % **************************************
    % YAHOO! INTEGER SOLUTION FOUND
    % **************************************
    
    if isempty(non_integer_binary) & isempty(non_integer_integer)  & isempty(non_semivar_semivar) & ~(output.problem == -1)
        if (cost<upper) & feasible
            x_min = x;
            upper = cost;
            [stack,lower] = prune(stack,upper,p.options,solved_nodes,p);
            
            if isfield(p.options,'plottruss')
                if p.options.plottruss
                    plottruss(3,'Best binary solution',poriginal,x_min);
                end
            end
            
        end       
        p = adaptivestrategy(p,upper,solved_nodes);        
        keep_digging = 0;
    end
    
    % **************************************
    % Stop digging if it won't give sufficient improvement anyway
    % **************************************
    if cost>upper*(1-1e-6)
        keep_digging = 0;
    end
    
    % **********************************
    % CONTINUE SPLITTING?
    % **********************************
    if keep_digging & (cost<upper)
        
        if solved_nodes == 1
            RootNodeInfeas =  TotalIntegerInfeas+TotalBinaryInfeas;
            RootNodeCost = cost;
        end
                
        % **********************************
        % BRANCH VARIABLE
        % **********************************
        [index,whatsplit,globalindex] = branchvariable(x,integer_variables,binary_variables,p.options,x_min,[],p);
        
        % **********************************
        % CREATE NEW PROBLEMS
        % **********************************
        p0_feasible = 1;
        p1_feasible = 1;
                
        switch whatsplit
            case 'binary'
                [p0,p1,index] = binarysplit(p,x,index,cost,[],sosgroups,sosvariables);
                
            case 'integer'
                [p0,p1] = integersplit(p,x,index,cost,x_min);
                
            case 'semi'
                [p0,p1] = semisplit(p,x,index,cost,x_min);
                
            case 'sos1'
                [p0,p1] = sos1split(p,x,index,cost,x_min);
                
            otherwise
        end
        
        node1.lb = p1.lb;
        node1.ub = p1.ub;
        node1.depth = p1.depth;
        node1.lower = p1.lower;
        node1.fixedvariable = globalindex;
        node1.fixdir = 'up';
        node1.TotalIntegerInfeas = TotalIntegerInfeas;
        node1.TotalBinaryInfeas = TotalBinaryInfeas;
        node1.IntInfeas = 1-(x(globalindex)-floor(x(globalindex)));
        node1.x0 = p1.x0;
        node1.binary_variables = p1.binary_variables;
        node1.semicont_variables = p1.semicont_variables;
        node1.semibounds = p1.semibounds;
        node1.pid = pid;pid = pid + 1;
        node1.sosgroups = p1.sosgroups;
        node1.sosvariables = p1.sosvariables;
        
        node0.lb = p0.lb;
        node0.ub = p0.ub;
        node0.depth = p0.depth;
        node0.lower = p0.lower;
        node0.fixedvariable = index;
        node0.fixdir = 'down';
        node0.TotalIntegerInfeas = TotalIntegerInfeas;
        node0.TotalBinaryInfeas = TotalBinaryInfeas;
        node0.IntInfeas = x(globalindex)-floor(x(globalindex));
        node0.x0 = p0.x0;
        node0.binary_variables = p0.binary_variables;
        node0.semicont_variables = p0.semicont_variables;
        node0.semibounds = p0.semibounds;
        node0.pid = pid;pid = pid + 1;
        node0.sosgroups = p0.sosgroups;
        node0.sosvariables = p0.sosvariables;
           
        if p0_feasible
            stack = push(stack,node0);
        end        
        if p1_feasible
            stack = push(stack,node1);
        end

    end
    
    % Lowest cost in any open node
    if ~isempty(stack)
        lower = min([stack.lower]);
        if can_use_ceil_lower
            lower = ceil(lower);
        end
    end
    
    % **********************************
    % Get a new node to solve
    % **********************************
    [node,stack] = pull(stack,p.options.bnb.method,x_min,upper);
    if ~isempty(node)
        p.lb = node.lb;
        p.ub = node.ub;
        p.depth = node.depth;
        p.lower = node.lower;
        p.fixedvariable = node.fixedvariable;
        p.fixdir = node.fixdir;
        p.TotalIntegerInfeas = node.TotalIntegerInfeas;
        p.TotalBinaryInfeas = node.TotalBinaryInfeas;
        p.IntInfeas =  node.IntInfeas;
        p.x0 = node.x0;
        p.binary_variables = node.binary_variables;
        p.semicont_variables = node.semicont_variables;
        p.semibounds = node.semibounds;
        % p.Musts = node.Musts;
        p.pid = node.pid;
        p.sosgroups = node.sosgroups;
        p.sosvariables = node.sosvariables;
    end
    gap = abs((upper-lower)/(1e-3+abs(upper)+abs(lower)));
    if isnan(gap)
        gap = inf;
    end
    
    if p.options.bnb.plotbounds
        lowerhist = [lowerhist lower];
        upperhist = [upperhist upper];
        hold off
        plot([lowerhist' upperhist']);
        drawnow
    end
    
    %DEBUG    if p.options.bnb.verbose;fprintf(' %4.0f : %12.3E  %7.2f   %12.3E  %2.0f   %2.0f %2.0f %2.0f %2.0f\n',solved_nodes,upper,100*gap,lower,length(stack)+length(node),sedd);end
    if p.options.bnb.verbose;fprintf(' %4.0f : %12.3E  %7.2f   %12.3E  %2.0f  %s\n',solved_nodes,upper,100*gap,lower,length(stack)+length(node),yalmiperror(output.problem));end
    
end
if p.options.bnb.verbose;showprogress([num2str2(solved_nodes,3)  ' Finishing.  Cost: ' num2str(upper) ],p.options.bnb.verbose);end


function stack = push(stackin,p)
if ~isempty(stackin)
    stack = [p;stackin];
else
    stack(1)=p;
end

%%
function [p,stack] = pull(stack,method,x_min,upper);

if ~isempty(stack)
    switch method
        case {'depth','depthfirst','depthbreadth','depthproject','depthbest'}
            [i,j]=max([stack.depth]);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);
            
        case 'project'
            [i,j]=min([stack.projection]);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);
            
        case 'breadth'
            [i,j]=min([stack.depth]);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);
            
        case 'best'
            [i,j]=min([stack.lower]);
            % candidates = find([stack.lower] == stack(j).lower);
            % [i,j] = min([stack(candidates).IntInfeas]);
            % j = candidates(j);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);
            
        otherwise
    end
else
    p = [];
end

% **********************************
%% BRANCH VARIABLE
% **********************************
function [index,whatsplit,globalindex] = branchvariable(x,integer_variables,binary_variables,options,x_min,Weight,p)
all_variables = [integer_variables(:);binary_variables(:)];

if isempty(setdiff(all_variables,p.sosvariables)) & strcmp(options.bnb.branchrule,'sos')
    % All variables are in SOS1 constraints
    for i = 1:length(p.sosgroups)        
        dist(i) = (sum(x(p.sosgroups{i}))-max(x(p.sosgroups{i})))/length(p.sosgroups{i});
    end
    % Which SOS to branch on
    [val,index] = max(dist);    
    whatsplit = 'sos1';
    globalindex = index;
    return
end

switch options.bnb.branchrule
    case 'weight'
        interror = abs(x(all_variables)-round(x(all_variables)));
        [val,index] = max(abs(p.weight(all_variables)).*interror);
    case 'first'
        index = min(find(abs(x(all_variables)-round(x(all_variables)))>options.bnb.inttol));
    case 'last'
        index = max(find(abs(x(all_variables)-round(x(all_variables)))>options.bnb.inttol));
    case 'min'
        nint = find(abs(x(all_variables)-round(x(all_variables)))>options.bnb.inttol);
        [val,index] = min(abs(x(nint)));
        index = nint(index);
    case 'max'
        [val,index] = max((abs(x(all_variables)-round(x(all_variables)))));
        %[val,index] = max(abs(p.c(all_variables)).^2.*(abs(x(all_variables)-round(x(all_variables)))));
    otherwise
        error('Branch-rule not supported')
end
if index<=length(integer_variables)
    whatsplit = 'integer';
    globalindex = integer_variables(index);
else
    index = index-length(integer_variables);
    whatsplit = 'binary';
    globalindex = binary_variables(index);
end
if isempty(index) | ~isempty(p.semicont_variables)
    for i = 1:length(p.semicont_variables)
        j = p.semicont_variables(i);
        if x(j)>= p.semibounds.lb(i) & x(j)<= p.semibounds.ub(i)
            s(i) = 0;
        elseif x(j)==0
            s(i) = 0;
        else
            s(i) = min([abs(x(j)-0); abs(x(j)-p.semibounds.lb(i));abs(x(j)-p.semibounds.ub(i))]);
        end
    end
    [val2,index2] = max(s);
    if isempty(val)
        whatsplit = 'semi';
        index = index2;
    elseif val2>val
        % index = p.semicont_variables(index);
        whatsplit = 'semi';
        index = index2;
    end
end



% **********************************
% SPLIT PROBLEM
% **********************************
function [p0,p1,variable] = binarysplit(p,x,index,lower,options,sosgroups,sosvariables)
p0 = p;
p1 = p;

variable = p.binary_variables(index);
tf = ~(ismembcYALMIP(p0.binary_variables,variable));
new_binary = p0.binary_variables(tf);

friends = [];
if ~isempty(sosvariables)
    if ismember(variable,sosvariables)
        i = 1;
        while i<=length(sosgroups)
            
            if ismember(variable,sosgroups{i})
                friends = setdiff(sosgroups{i},variable);
                break
            else
                i = i + 1;
            end
        end
    end
end

p0.ub(variable)=0;
p0.lb(variable)=0;
if length(friends) == 1
     p0.ub(friends) = 1;
     p0.lb(friends) = 1;
end

p0.lower = lower;
p0.depth = p.depth+1;
p0.binary_variables = new_binary;%setdiff1D(p0.binary_variables,variable);
%p0.binary_variables = setdiff(p0.binary_variables,friends);

p1.ub(variable)=1;
p1.lb(variable)=1;
if length(friends) > 1
    p1.ub(friends)=0;
    p1.lb(friends)=0;
end

p1.binary_variables = new_binary;%p0.binary_variables;%setdiff1D(p1.binary_variables,variable);
%p1.binary_variables = setdiff(p1.binary_variables,friends);
p1.lower = lower;
p1.depth = p.depth+1;

% % *****************************
% % PROCESS MOST PROMISING FIRST
% % (p0 in top of stack)
% % *****************************
if x(variable)>0.5
    pt=p1;
    p1=p0;
    p0=pt;
end

function [p0,p1] = integersplit(p,x,index,lower,options,x_min)

variable = p.integer_variables(index);
current = x(p.integer_variables(index));
lb = floor(current)+1;
ub = floor(current);

% xi<ub
p0 = p;
p0.lower = lower;
p0.depth = p.depth+1;
p0.x0(variable) = ub;
p0.ub(variable)=min(p0.ub(variable),ub);

% xi>lb
p1 = p;
p1.lower = lower;
p1.depth = p.depth+1;
p1.x0(variable) = lb;
p1.lb(variable)=max(p1.lb(variable),lb);

% *****************************
% PROCESS MOST PROMISING FIRST
% *****************************
if lb-current<0.5
    pt=p1;
    p1=p0;
    p0=pt;
end

function [p0,p1] = sos1split(p,x,index,lower,options,x_min)

v = p.sosgroups{index};
n = ceil(length(v)/2);
v1 = v(randperm(length(v),n));
v2 = setdiff(v,v1);
%v1 = v(1:n);
%v2 = v(n+1:end);

% In first node, set v2 to 0 and v1 to sosgroup
p0 = p;p0.lower = lower;
p0.sosgroups{index} = v1;
p0.ub(v2) = 0;

% In second node, set v1 to 0 and v1 to sosgroup
p1 = p;p1.lower = lower;
p1.sosgroups{index} = v2;
p1.ub(v1) = 0;


function [p0,p1] = semisplit(p,x,index,lower,options,x_min)

variable = p.semicont_variables(index);
current = x(p.semicont_variables(index));

p0 = p;
p0.lower = lower;
p0.depth = p.depth+1;
p0.x0(variable) = 0;
p0.lb(variable)=0;
p0.ub(variable)=0;

p1 = p;
p1.lower = lower;
p1.depth = p.depth+1;
p1.x0(variable) = p.semibounds.lb(index);
p1.lb(variable) = p.semibounds.lb(index);
p1.ub(variable) = p.semibounds.ub(index);

p0.semicont_variables = setdiff(p.semicont_variables,variable);
p1.semicont_variables = setdiff(p.semicont_variables,variable);
p0.semibounds.lb(index)=[];
p0.semibounds.ub(index)=[];
p1.semibounds.lb(index)=[];
p1.semibounds.ub(index)=[];



function s = num2str2(x,d,c);
if nargin==3
    s = num2str(x,c);
else
    s = num2str(x);
end
s = [repmat(' ',1,d-length(s)) s];


function [stack,lower] = prune(stack,upper,options,solved_nodes,p)
% *********************************
% PRUNE STACK W.R.T NEW UPPER BOUND
% *********************************
if ~isempty(stack)
    %    toolarge = find([stack.lower]>upper*(1-1e-4));
    toolarge = find([stack.lower]>upper*(1-options.bnb.prunetol));
    if ~isempty(toolarge)
        stack(toolarge)=[];
    end
end

if ~isempty(stack)
    lower = min([stack.lower]);
else
    lower = upper;
end

function p = adaptivestrategy(p,upper,solved_nodes)
% **********************************'
% SWITCH NODE SELECTION STRATEGY?
% **********************************'
if strcmp(p.options.bnb.method,'depthproject') & (upper<inf)
    p.options.bnb.method = 'project';
end
if strcmp(p.options.bnb.method,'depthbest') & (upper<inf)
    p.options.bnb.method = 'best';
end
if strcmp(p.options.bnb.method,'depthprojection') & (upper<inf)
    p.options.bnb.method = 'projection';
end
if strcmp(p.options.bnb.method,'depthbreadth') & (upper<inf)
    p.options.bnb.method = 'breadth';
end
if strcmp(p.options.bnb.method,'depthest') & (upper<inf)
    p.options.bnb.method = 'est';
end

function res = resids(p,x)
res= [];
if p.K.f>0
    res = -abs(p.F_struc(1:p.K.f,:)*[1;x]);
end
if p.K.l>0
    res = [res;p.F_struc(p.K.f+1:p.K.f+p.K.l,:)*[1;x]];
end
if (length(p.K.s)>1) | p.K.s>0
    top = 1+p.K.f+p.K.l;
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        X = p.F_struc(top:top+n^2-1,:)*[1;x];top = top+n^2;
        X = reshape(X,n,n);
        res = [res;min(eig(X))];
    end
end
res = [res;min([p.ub-x;x-p.lb])];

function p = Updatecostbound(p,upper,lower);
if p.simplecost
    if ~isinf(upper)
        ind = find(p.c);
        if p.c(ind)>0
            p.ub(ind) = min(p.ub(ind),(upper-p.f)/p.c(ind));
        else
            p.lb(ind) = max(p.lb(ind),(p.f-upper)/abs(p.c(ind)));
        end
    end
end

function [x_min,upper] = initializesolution(p);

x_min = zeros(length(p.c),1);
upper = inf;
if p.options.usex0
    z = p.x0;
    residual = resids(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-1e-12) & all(residual(1+p.K.f:end)>=-1e-6);
    if relaxed_feasible & all(z(p.integer_variables)==fix(z(p.integer_variables))) & all(z(p.binary_variables)==fix(z(p.binary_variables)))
        upper = computecost(p.f,p.corig,p.Q,z,p);%upper = p.f+p.c'*z+z'*p.Q*z;
        x_min = z;
    end
else
    p.x0 = zeros(length(p.c),1);
    x = p.x0;
    z = evaluate_nonlinear(p,x);
    residual = resids(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=p.options.bmibnb.pdtol);
    if relaxed_feasible
        upper = computecost(p.f,p.corig,p.Q,z,p);%upper = p.f+p.c'*z+z'*p.Q*z;
        x_min = x;
    end
end



function [p,poriginal,stack] = pruneglobally(p,poriginal,upper,lower,stack,x);

if isempty(p.nonlinear) & (nnz(p.Q)==0) & p.options.bnb.nodetight
    pp = poriginal;
    
    if p.K.l > 0
        A = -pp.F_struc(1+pp.K.f:pp.K.f+pp.K.l,2:end);
        b = pp.F_struc(1+p.K.f:p.K.f+p.K.l,1);
    else
        A = [];
        b = [];
    end
    
    if (nnz(p.Q)==0) & ~isinf(upper)
        A = [pp.c';-pp.c';A];
        b = [upper;-(lower-0.0001);b];
    else
        % c = p.c;
        % Q = p.Q;
        % A = [c'+2*x'*Q;A];
        % b = [2*x'*Q*x+c'*x;b];
    end
    
    [lb,ub,redundant,pss] = milppresolve(A,b,pp.lb,pp.ub,pp.integer_variables,pp.binary_variables,ones(length(pp.lb),1));
    
    if ~isempty(redundant)
        if (nnz(p.Q)==0) & ~isinf(upper)
            redundant = redundant(redundant>2)-2;
        else
            %    redundant = redundant(redundant>1)-1;
        end
        if length(redundant)>0
            poriginal.K.l=poriginal.K.l-length(redundant);
            poriginal.F_struc(poriginal.K.f+redundant,:)=[];
            p.K.l=p.K.l-length(redundant);
            p.F_struc(p.K.f+redundant,:)=[];
        end
    end
    if ~isempty(stack)
        keep = ones(length(stack),1);
        for i = 1:length(stack)
            stack(i).lb = max([stack(i).lb lb]')';
            stack(i).ub = min([stack(i).ub ub]')';
            if any(stack(i).lb>stack(i).ub)
                keep(i) = 0;
            end
        end
        stack = stack(find(keep));
    end
    poriginal.lb = max([poriginal.lb lb]')';
    poriginal.ub = min([poriginal.ub ub]')';
    p.lb = max([p.lb lb]')';
    p.ub = min([p.ub ub]')';
end


function [p,poriginal,stack] = fixvariables(p,poriginal,upper,lower,stack,x_min,monotinicity)
% Fix variables

if p.options.bnb.nodefix & (p.K.f == 0) & (nnz(p.Q)==0) & isempty(p.nonlinear)
    
    A = -poriginal.F_struc(poriginal.K.f + (1:poriginal.K.l),2:end);
    b = poriginal.F_struc(poriginal.K.f + (1:poriginal.K.l),1);
    c = poriginal.c;
    [fix_up,fix_down] = presolve_fixvariables(A,b,c,poriginal.lb,poriginal.ub,monotinicity);
    %
    poriginal.lb(fix_up) = 1;
    p.lb(fix_up) = 1;
    
    %     not_in_obj = find(p.c==0);
    %     constrained_blow = all(poriginal.F_struc(1:poriginal.K.l,1+not_in_obj)>=0,1);
    %     sdp_positive = sdpmonotinicity(not_in_obj) == -1;
    %     can_fix = not_in_obj(find(constrained_blow & sdp_positive));
    %
    %     still_on = find(p.lb==0 & p.ub==1);
    %     p.lb(intersect(can_fix,still_on)) = 1;
    %     still_on = find(poriginal.lb==0 & poriginal.ub==1);
    %     poriginal.lb(intersect(can_fix,still_on)) = 1;
    
    if ~isempty(stack) & ~(isempty(fix_up)  & isempty(fix_down))
        keep = ones(length(stack),1);
        for i = 1:length(stack)
            stack(i).lb = max([stack(i).lb poriginal.lb]')';
            stack(i).ub = min([stack(i).ub poriginal.ub]')';
            if any(stack(i).lb>stack(i).ub)
                keep(i) = 0;
            end
        end
        stack = stack(find(keep));
    end
end


function [feasible,p] = checkmusts(p)
feasible = 1;

if ~isempty(p.Musts)
    
    %     for i = 1:size(p.Musts,1)
    %         if all(p.ub(find(p.Musts(i,:)))==0)
    %             1
    %         end
    %     end
    %
    %     TurnedOff = find(p.ub == 0);
    %     if ~isempty(TurnedOff)
    %         p.Musts(:,TurnedOff) = 0;
    %         Failure = find(sum(p.Musts,2)==0);
    %         if ~isempty(Failure)
    %             1;%feasible = 0;
    %             return
    %         end
    %         TurnedOn = find(sum(p.Musts,2)==1);
    %         if ~isempty(TurnedOn)
    %             %    p.lb(TurnedOn) = 1;
    %         end
    %     end
end


function stack = prunecardinality(p,poriginal,stack,lower,upper)
% if length(poriginal.binary_variables)==length(poriginal.c) & nnz(poriginal.Q)==0 & all(poriginal.c)>0
%     card_max_low = lower/max(p.c);
%     card_max_high = upper/min(p.c);
%     [card_max_low card_max_high]
%     keep = ones(1,length(stack));
%     for i = 1:length(stack)
%         fixed_zero = nnz(stack(i).ub==0);
%         fixed_one  = nnz(stack(i).lb==1);
%         if fixed_one>= card_max_high
%             1
%         end
%     end
% end




