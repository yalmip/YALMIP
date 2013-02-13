function [sol,m,Q,residuals,everything] = solvesos(F,obj,options,params,candidateMonomials)
%SOLVESOS Sum of squares decomposition
%
%    [sol,m,B,residual] = solvesos(F,h,options,params,monomials) is used
%    for finding SOS decompositions of polynomials.
%
%    The coefficients of the polynomials are assumed linear w.r.t a set of
%    decision variables params and polynomial with respect to a variable x.
%    
%    An extension with a nonlinear parameterization in params is possible.
%    Note though that this gives BMIs or PMIs, solvable (locally) only if
%    PENBMI is installed, or by specifying 'moment' as solver to try to
%    solve the nonconvex semidefinite programming problem using a
%    semidefinite relaxation based on moments.
%
%    The SOS problem can be formulated as
%
%              min h(params)
%
%       subject to  F(i) >(=) 0 or F(i) is SOS w.r.t x
%
%   INPUT
%    F         : SET object with SOS constrained polynomials and constraints on variables params
%    h         : scalar SDPVAR object (can be [])
%    options   : options structure obtained from SDPSETTINGS (can be [])
%    params    : SDPVAR object defining parametric variables (can be [])
%    monomials : SDPVAR object with user-specified monomials for decomposition (can be [])
%
%   OUTPUT
%    sol       : Solution diagnostic from SDP problem
%    v         : Cell with monomials used in decompositions
%    Q         : Cell with Gram matrices, p = v{i}'*Q{i}*v{i}, where p is the ith SOS polynomial in your model.
%    residuals : Mismatch between p and decompositions. Same values (modulo numerical issue) as checkset(find(is(F,'sos')))
%                Warning, these residuals are not computed on matrix sos-of-squares
%
%   EXAMPLE
%    x = sdpvar(1);solvesos(set(sos(x^4+x^3+1)));                    % Simple decompositions
%    x = sdpvar(1);t = sdpvar(1);solvesos(set(sos(x^4+x^3+1-t)),-t); % Lower bound by maximizing t
%
%   NOTES
%
%    Variables not part of params, but part of non-SOS constraints in F
%    or objective h will automatically be appended to the params list.
%
%    To extract SOS decomposition, use the command SOSD (or compute from use v and Q)
%
%    If the 5th input argument is used, no additional monomial reduction is
%    performed (Newton, inconstency, congruence). It is thus assumed that
%    the supplied candidate monomials constitute a sufficient basis.
%
%    The field options.sos can be used to tune the SOS-calculations. See HTML help for details
%
%     sos.model          - Kernel (1) or image (2) representation of SOS problem [0|1|2 (0, YALMIP decides)]
%     sos.newton         - Use Newton polytope to reduce size [0|1 (1)]
%     sos.congruence     - Block-diagonalize using congruence classes [0|1|2 (2)]
%     sos.scale          - Scale polynomial [0|1 (1)]
%     sos.numblkdg       - Try to perform a-posteriori block-diagonalization [real  (0)]
%     sos.inconsistent   - Remove diagonal-inconsistent monomials [0|1|2 (0)]
%     sos.clean          - Remove monomials with coefficients < clean [real > 0 (1e-4)]
%     sos.traceobj       - Minimize trace of Gram matrix in problems without objective function [0|1 (0)]
%     sos.extlp          - Extract simple translated LP cones when performing dualization [0|1 (1)]
%
% See also SOS, SOSD, SDPSETTINGS, SOLVEMOMENT, SDPVAR, SDISPLAY

%% Time YALMIP
yalmip_time = clock;

% ************************************************
%% Check #inputs
% ************************************************
if nargin<5
    candidateMonomials = [];
    if nargin<4
        params = [];
        if nargin<3
            options = sdpsettings;
            if nargin<2
                obj = [];
                if nargin<1
                    help solvesos
                    return
                end
            end
        end
    end
end

if isa(obj,'double')
    obj = [];
end

if isa(F,'constraint')
    F = set(F);
end

if isequal(obj,0)
    obj = [];
end

if isempty(options)
    options = sdpsettings;
end

% Lazy syntax (not official...)
if nargin==1 & isa(F,'sdpvar')
    F = set(sos(F));
end

if ~isempty(options)
    if options.sos.numblkdg
        [sol,m,Q,residuals,everything] = solvesos_find_blocks(F,obj,options,params,candidateMonomials);
        return
    end
end

% *************************************************************************
%% Extract all SOS constraints and candidate monomials
% *************************************************************************
if ~any(is(F,'sos'))
    error('At-least one constraint should be an SOS constraints!');
end
p = [];
ranks = [];
for i = 1:length(F)
    if is(F(i),'sos')
        pi = sdpvar(F(i));
        p{end+1} = pi;
        ranks(end+1) = getsosrank(pi); % Desired rank : Experimental code
    end
end
if isempty(candidateMonomials)
    for i = 1:length(F)
        candidateMonomials{i}=[];
    end
elseif isa(candidateMonomials,'sdpvar')
    cM=candidateMonomials;
    candidateMonomials={};
    for i = 1:length(p)
        candidateMonomials{i}=cM;
    end
elseif isa(candidateMonomials,'cell')
    if length(p)~=length(candidateMonomials)
        error('Dimension mismatch between the candidate monomials and the number of SOS constraints');
    end
end

% *************************************************************************
%% Get the parametric constraints
% *************************************************************************
F_original = F;
F_parametric = F(find(~is(F,'sos')));
if isempty(F_parametric)
    F_parametric = set([]);
end
% Expand removes so called marker variables generated in ismember
ParametricBefore = getvariables(F_parametric);

% *************************************************************************
%% Expand the parametric constraints
% *************************************************************************
if ~isempty(yalmip('extvariables'))
    [F_parametric,failure] = expandmodel(F_parametric,obj,options);
    F_parametric = expanded(F_parametric,1);
    obj = expanded(obj,1);    
    if failure
        Q{1} = [];m{1} = [];residuals = [];everything = [];
        sol.yalmiptime = etime(clock,yalmip_time);
        sol.solvertime = 0;
        sol.info = yalmiperror(14,'YALMIP');
        sol.problem = 14;
    end
end

if ~isempty(params)
    if ~isa(params,'sdpvar')
        error('Fourth argment should be a SDPVAR variable or empty')
    end
end

% *************************************************************************
% Collect all possible parametric variables
% *************************************************************************
ParametricVariables = uniquestripped([depends(obj) depends(F_parametric) depends(params) ParametricBefore]);

if any(find(is(F_parametric,'parametric')))
    F_parametric(find(is(F_parametric,'parametric')))=[];
end
if any(find(is(F,'parametric')))
    F(find(is(F,'parametric')))=[];
end

if options.verbose>0;
    disp('-------------------------------------------------------------------------');
    disp('YALMIP SOS module started...');
    disp('-------------------------------------------------------------------------');
end

% *************************************************************************
%% INITIALIZE SOS-DECOMPOSITIONS SDP CONSTRAINTS
% *************************************************************************
F_sos = set([]);

% *************************************************************************
%% FIGURE OUT ALL USED PARAMETRIC VARIABLES
% *************************************************************************
AllVariables =  uniquestripped([depends(obj) depends(F_original) depends(F_parametric)]);
ParametricVariables = intersect(ParametricVariables,AllVariables);
MonomVariables = setdiff(AllVariables,ParametricVariables);
params = recover(ParametricVariables);
if isempty(MonomVariables)
    error('No independent variables? Perhaps you added a constraint set(p(x)) when you meant set(sos(p(x)))');
end
if options.verbose>0;disp(['Detected ' num2str(length(ParametricVariables)) ' parametric variables and ' num2str(length(MonomVariables)) ' independent variables.']);end

% ************************************************
%% ANY BMI STUFF
% ************************************************
NonLinearParameterization = 0;
if ~isempty(ParametricVariables)
    monomtable = yalmip('monomtable');
    ParametricMonomials = monomtable(uniquestripped([getvariables(obj) getvariables(F_original)]),ParametricVariables);
    if any(sum(abs(ParametricMonomials),2)>1)
        NonLinearParameterization = 1;
    end
end

% ************************************************
%% ANY INTEGER DATA
% ************************************************
IntegerData = 0;
if ~isempty(ParametricVariables)
    globalInteger =  [yalmip('binvariables') yalmip('intvariables')];    
    integerVariables = getvariables(F_parametric(find(is(F_parametric,'binary') | is(F_parametric,'integer'))));
    integerVariables = [integerVariables intersect(ParametricVariables,globalInteger)];
    integerVariables = intersect(integerVariables,ParametricVariables);
    IntegerData = ~isempty(integerVariables);
end

% ************************************************
%% ANY UNCERTAIN DATA
% ************************************************
UncertainData = 0;
if ~isempty(ParametricVariables)
    UncertainData = any(is(F_parametric,'uncertain'));  
end

% ************************************************
%% Any Interval-data
% ************************************************
IntervalData = any(is(F,'interval')) |  any(is(F_parametric,'interval'));
if ~isempty(obj)
     IntervalData | is(obj,'interval');
end

% ************************************************
%% DISPLAY WHAT WE FOUND
% ************************************************
if options.verbose>0 & ~isempty(F_parametric)
    nLP = 0;
    nEQ = 0;
    nLMI = sum(full(is(F_parametric,'lmi')) &  full(~is(F_parametric,'element-wise'))); %FULL due to bug in ML 7.0.1
    for i = 1:length(F_parametric)
        if is(F_parametric,'element-wise')
            nLP = nLP + prod(size(F_parametric(i)));
        end
        if is(F_parametric,'equality')
            nEQ = nEQ + prod(size(F_parametric(i)));
        end
    end
    disp(['Detected ' num2str(full(nLP)) ' linear inequalities, ' num2str(full(nEQ)) ' equality constraints and ' num2str(full(nLMI)) ' LMIs.']);
end

% ************************************************
%% IMAGE OR KERNEL REPRESENTATION?
% ************************************************
noRANK = all(isinf(ranks));
switch options.sos.model
    case 0
        constraint_classes = constraintclass(F);
        noCOMPLICATING = ~any(ismember([7 8 9 10 12 14 15],constraint_classes));
        if noCOMPLICATING & ~NonLinearParameterization & noRANK & ~IntegerData & ~IntervalData
            options.sos.model = 1;
            if options.verbose>0;disp('Using kernel representation (options.sos.model=1).');end
        else
            if NonLinearParameterization
                if options.verbose>0;disp('Using image representation (options.sos.model=2). Nonlinear parameterization found');end
            elseif ~noRANK
                if options.verbose>0;disp('Using image representation (options.sos.model=2). SOS-rank constraint was found.');end
            elseif IntegerData
                if options.verbose>0;disp('Using image representation (options.sos.model=2). Integrality constraint was found.');end
            elseif UncertainData
                if options.verbose>0;disp('Using image representation (options.sos.model=2). Uncertain data was found.');end                
            elseif IntervalData
                if options.verbose>0;disp('Using image representation (options.sos.model=2). Interval-data was found.');end                                
            else
                if options.verbose>0;disp('Using image representation (options.sos.model=2). Integer data, KYPs or similar was found.');end
            end
            options.sos.model = 2;
        end
    case 1
        if NonLinearParameterization
            if options.verbose>0;disp('Switching to image model due to nonlinear parameterization (not supported in kernel model).');end
            options.sos.model = 2;
        end
        if ~noRANK
            if options.verbose>0;disp('Switching to image model due to SOS-rank constraints (not supported in kernel model).');end
            options.sos.model = 2;
        end
        if IntegerData
            if options.verbose>0;disp('Switching to image model due to integrality constraints (not supported in kernel model).');end
            options.sos.model = 2;
        end        
    case 3
    otherwise
end

if ~isempty(yalmip('extvariables')) & options.sos.model == 2 & nargin<4
    disp(' ')
    disp('**Using nonlinear operators in SOS problems can cause problems.')
    disp('**Please specify all parametric variables using the fourth argument');
    disp(' ');
end

% ************************************************
%% SKIP DIAGONAL INCONSISTENCY FOR PARAMETRIC MODEL
% ************************************************
if ~isempty(params) & options.sos.inconsistent
    if options.verbose>0;disp('Turning off inconsistency based reduction (not supported in parametric models).');end
    options.sos.inconsistent = 0;
end

% ************************************************
%% INITIALIZE OBJECTIVE
% ************************************************
if ~isempty(obj)
    options.sos.traceobj = 0;
end
parobj = obj;
obj    = [];

% ************************************************
%% SCALE SOS CONSTRAINTS 
% ************************************************
if options.sos.scale
    for constraint = 1:length(p)
        normp(constraint) = sqrt(norm(full(getbase(p{constraint}))));
        p{constraint} = p{constraint}/normp(constraint);
        sizep(constraint) = size(p{constraint},1);
    end
else
    normp = ones(length(p),1);
end

% ************************************************
%% Some stuff not supported for
%  matrix valued SOS yet, turn off for safety
% ************************************************
for constraint = 1:length(p)
    sizep(constraint) = size(p{constraint},1);
end
if any(sizep>1)
    options.sos.postprocess = 0;
    options.sos.reuse = 0;
end

% ************************************************
%% SKIP CONGRUENCE REDUCTION WHEN SOS-RANK
% ************************************************
if ~all(isinf(ranks))
    options.sos.congruence = 0;
end

% ************************************************
%% Create an LP model to speed up things in Newton
%  polytope reduction
% ************************************************
if options.sos.newton
    temp=sdpvar(1,1);
    tempops = options;
    tempops.solver = 'cdd,glpk,*';  % CDD is generally robust on these problems
    tempops.verbose = 0;
    tempops.saveduals = 0;
    [aux1,aux2,aux3,LPmodel] = export(set(temp>=0),temp,tempops);   
else
    LPmodel = [];
end


% ************************************************
%% LOOP THROUGH ALL SOS CONSTRAINTS
% ************************************************
for constraint = 1:length(p)
    % *********************************************
    %% FIND THE VARIABLES IN p, SORT, UNIQUE ETC
    % *********************************************
    if options.verbose>1;disp(['Creating SOS-description ' num2str(constraint) '/' num2str(length(p)) ]);end

    pVariables = depends(p{constraint});
    AllVariables = uniquestripped([pVariables ParametricVariables]);
    MonomVariables = setdiff1D(pVariables,ParametricVariables);
    x = recover(MonomVariables);
    z = recover(AllVariables);
    MonomIndicies = find(ismember(AllVariables,MonomVariables));
    ParametricIndicies = find(ismember(AllVariables,ParametricVariables));

    if isempty(MonomIndicies)
        % This is the case set(sos(t)) where t is a parametric (matrix) variable
        % This used to create an error message befgore to avoid some silly
        % bug in the model generation. Creating this error message is
        % stupid, but at the same time I can not remember where the bug was
        % and I have no regression test  for this case. To avoid
        % introducing same bug again by mistake, I create all data
        % specifically for this case        
        previous_exponent_p_monoms = [];%exponent_p_monoms;
        n = length(p{constraint});
        A_basis = getbase(sdpvar(n,n,'full'));d = find(triu(ones(n)));A_basis = A_basis(d,2:end);
        BlockedA{constraint} = {A_basis};
        Blockedb{constraint} = p{constraint}(d);        
        BlockedN{constraint} = {zeros(1,0)};
        Blockedx{constraint} = x;
        Blockedvarchange{constraint}=zeros(1,0);    
        continue
        %  error('You have constraints of the type set(sos(f(parametric_variables))). Please use set(f(parametric_variables) > 0) instead')
    end

    % *********************************************
    %% Express p in monimials and coefficients  
    % *********************************************
    [exponent_p,p_base] = getexponentbase(p{constraint},z);
    
    % *********************************************
    %% Powers for user defined candidate monomials
    % (still experimental)
    % *********************************************
    if ~all(cellfun('isempty',candidateMonomials))
        exponent_c = [];
        if isa(candidateMonomials{constraint},'cell')            
            for i = 1:length(candidateMonomials{constraint})
                 exponent_c{i} = getexponentbase(candidateMonomials{constraint}{i},z);
                 exponent_c{i} = exponent_c{i}(:,MonomIndicies);
            end
        else
            exponent_c{1} = getexponentbase(candidateMonomials{constraint},z);
            exponent_c{1} = exponent_c{1}(:,MonomIndicies);
        end
    else
        exponent_c = [];
    end
       
    % *********************************************
    %% STUPID PROBLEM WITH ODD HIGHEST POWER?...
    % *********************************************
    if isempty(ParametricIndicies)
        max_degrees = max(exponent_p(:,MonomIndicies),[],1);
        bad_max = any(max_degrees-fix((max_degrees/2))*2);
        if bad_max
            for i = 1:length(p)
                Q{i}=[];
                m{i}=[];
            end
            residuals=[];
            everything = [];
            sol.yalmiptime = etime(clock,yalmip_time);
            sol.solvertime = 0;
            sol.info = yalmiperror(1,'YALMIP');
            sol.problem = 2;
            return
        end
    end

    % *********************************************
    %% Can we make a smart variable change (no code)
    % *********************************************
    exponent_p_monoms = exponent_p(:,MonomIndicies);
    varchange = ones(1,size(MonomIndicies,2));

    % *********************************************
    %% Unique monoms (copies due to parametric terms)
    % *********************************************
    exponent_p_monoms = uniquesafe(exponent_p_monoms,'rows');

    if options.sos.reuse & constraint > 1 & isequal(previous_exponent_p_monoms,exponent_p_monoms)
        % We don't have to do anything, candidate monomials can be-used
        if options.verbose>1;disp(['Re-using all candidate monomials (same problem structure)']);end
    else       

        % *********************************************
        % User has supplied the whole candidate structure
        % Don't process this
        % *********************************************
        if ~isempty(exponent_c)
            exponent_m{1} = [];
            N = {};
            for i = 1:length(exponent_c)
                exponent_m{i} = [exponent_m{1};exponent_c{i}];
                N{i,1} = exponent_c{i};
            end
        else
            % *********************************************
            %% CORRELATIVE SPARSITY PATTERN
            % *********************************************
            [C,csclasses] = corrsparsity(exponent_p_monoms,options);

            % *********************************************
            %% GENERATE MONOMIALS
            % *********************************************
            exponent_m = monomialgeneration(exponent_p_monoms,csclasses);

            % *********************************************
            %% REDUCE #of MONOMIALS
            % *********************************************
            % Fix for matrix case, perform newton w.r.t
            % diagonal polynomials only. This can be
            % improved, but for now, keep it simple...
            n = length(p{constraint});diag_elements = 1:(n+1):n^2;used_diagonal = find(any(p_base(diag_elements,:),1));
            exponent_p_monoms_diag = exponent_p(used_diagonal,MonomIndicies);
            exponent_m = monomialreduction(exponent_m,exponent_p_monoms_diag,options,csclasses,LPmodel);

            % *********************************************
            %% BLOCK PARTITION THE MONOMIALS BY CONGRUENCE
            % *********************************************
            N = congruenceblocks(exponent_m,exponent_p_monoms,options,csclasses);
            % *********************************************
            %% REDUCE FURTHER BY EXPLOITING BLOCK-STRUCTURE
            % *********************************************
            N = blockmonomialreduction(exponent_p_monoms_diag,N,options);
            
        end


        % *********************************************
        %% PREPARE FOR SDP FORMULATION BY CALCULATING ALL
        % POSSIBLE MONOMIAL PRODUCS IN EACH BLOCK
        % *********************************************
        [exponent_m2,N_unique] = monomialproducts(N);

        % *********************************************
        %% CHECK FOR BUG/IDIOT PROBLEMS IN FIXED PROBLEM
        % *********************************************
        if isempty(ParametricIndicies)
            if ~isempty(setdiff(exponent_p_monoms,N_unique(:,3:end),'rows'))
                for i = 1:length(p)
                    Q{i} = [];
                    m{i} = [];
                end
                residuals = [];everything = [];
                sol.problem = 2;
                sol.info = yalmiperror(1,'YALMIP');
                warning('Problem is trivially infeasible (odd highest power?)');
                return
            end
        end
    end

    previous_exponent_p_monoms = exponent_p_monoms;

    % *********************************************
    %% GENERATE DATA FOR SDP FORMULATIONS
    % *********************************************
    p_base_parametric = [];
    n = length(p{constraint});
    for i=1:length(p{constraint})
        for j = 1:length(p{constraint})
            p_base_parametric = [p_base_parametric parameterizedbase(p{constraint}(i,j),z,params,ParametricIndicies,exponent_p,p_base((i-1)*n+j,:))];
        end
    end
    [BlockedA{constraint},Blockedb{constraint}] = generate_kernel_representation_data(N,N_unique,exponent_m2,exponent_p,p{constraint},options,p_base_parametric,ParametricIndicies,MonomIndicies,constraint==1);

    % SAVE FOR LATER
    BlockedN{constraint} = N;
    Blockedx{constraint} = x;
    Blockedvarchange{constraint}=varchange;
end

% *********************************************
%% And now get the SDP formulations
%
% The code above has generated matrices A and b
% in AQ == b(parametric)
%
% We use these to generate kernel or image models
% *********************************************
sol.problem = 0;
switch options.sos.model
    case 1
        % Kernel model
        [F,obj,BlockedQ,Primal_matrices,Free_variables] = create_kernelmodel(BlockedA,Blockedb,F_parametric,parobj,options,[]);
    case 2
        % Image model
        [F,obj,BlockedQ,sol] = create_imagemodel(BlockedA,Blockedb,F_parametric,parobj,options);

    case 3
        % Un-official model to solve bilinearly parameterized SOS using SDPLR
        [F,obj,options] = create_lrmodel(BlockedA,Blockedb,F_parametric,parobj,options,ParametricVariables);

    otherwise
end

% Unofficial fifth output with pseudo-numerical data
everything.BlockedA = BlockedA;
everything.Blockedb = Blockedb;
everything.F = F;
everything.obj = obj;

% % **********************************************
% %% SOLVE SDP
% % **********************************************
if sol.problem == 0
    if options.verbose > 0
        disp(' ');
    end
    if all(isinf(ranks))
        % Standard case
        %sol = solvesdp(F+set(-10<recover(depends(F)) < 10),obj,sdpsettings('bmibnb.maxit',200,'solver','bmibnb','bmibnb.lpred',1))
        sol =  solvesdp(F,obj,options);
    else
        % We have to alter the problem slightly if there are rank
        % constraints on the decompositions
        sol =  solveranksos(F,obj,options,ranks,BlockedQ);      
    end
end

% **********************************************
%% Recover solution (and rescale model+solution)
% **********************************************
for constraint = 1:length(p)
    for i = 1:length(BlockedQ{constraint})
        doubleQ{constraint}{i} = normp(constraint)*double(BlockedQ{constraint}{i});
    end
    doubleb{constraint}=normp(constraint)*double(Blockedb{constraint});
end

% **********************************************
%% Minor post-process
% **********************************************
if all(sizep<=1)
    [doubleQ,residuals] = postprocesssos(BlockedA,doubleb,doubleQ,[],options);
else
    residuals = 0;
end

% **********************************************
%% Safety check for bad solvers...
% **********************************************
if any(residuals > 1e-3) & (sol.problem == 0) & options.verbose>0
    disp(' ');
    disp('-> Although the solver indicates no problems,')
    disp('-> the residuals in the problem are really bad.')
    disp('-> My guess: the problem is probably infeasible.')
    disp('-> Make sure to check how well your decomposition')
    disp('-> matches your polynomial (see manual)')
    disp('-> You can also try to change the option sos.model')
    disp('-> or use another SDP solver.')
    disp(' ');
end

% **********************************************
%% Confused users. Primal dual kernel image?...
% **********************************************
if options.verbose > 0
    if sol.problem == 1
        if options.sos.model == 1
            disp(' ')
            disp('-> Solver reported infeasible dual problem.')
            disp('-> Your SOS problem is probably unbounded.')
        elseif options.sos.model == 2
            disp(' ')
            disp('-> Solver reported infeasible primal problem.')
            disp('-> Your SOS problem is probably infeasible.')
        end
    elseif sol.problem == 2
        if options.sos.model == 1
            disp(' ')
            disp('-> Solver reported unboundness of the dual problem.')
            disp('-> Your SOS problem is probably infeasible.')
        elseif options.sos.model == 2
            disp(' ')
            disp('-> Solver reported unboundness of the primal problem.')
            disp('-> Your SOS problem is probably unbounded.')            
        end
    elseif sol.problem == 12
            disp(' ')
            disp('-> Solver reported unboundness or infeasibility of the primal problem.')
            disp('-> Your SOS problem is probably unbounded.')            
    end
end

% **********************************************
%% De-block
% **********************************************
for constraint = 1:length(p)
    Qtemp = [];
    for i = 1:length(BlockedQ{constraint})
        Qtemp = blkdiag(Qtemp,doubleQ{constraint}{i});
    end
    Q{constraint} = full(Qtemp);
end

% **********************************************
%% Experimental code for declaring sparsity
% **********************************************
if options.sos.impsparse == 1
    somesmall = 0;
    for i = 1:length(BlockedQ)
        for j =  1:length(BlockedQ{i})
            small = find(abs(double(BlockedQ{i}{j}))<options.sos.sparsetol);
            nullThese{i}{j} = small;
            if ~isempty(small)
                somesmall = 1;
            end
        end
    end
    if somesmall
        [F,obj,BlockedQ,Primal_matrices,Free_variables] = create_kernelmodel(BlockedA,Blockedb,F_parametric,parobj,options,nullThese);
        sol =  solvesdp(F,obj,options);
        for constraint = 1:length(p)
            for i = 1:length(BlockedQ{constraint})
                doubleQ{constraint}{i} = normp(constraint)*double(BlockedQ{constraint}{i});
            end
            doubleb{constraint}=normp(constraint)*double(Blockedb{constraint});
        end

        % **********************************************
        %% Post-process
        % **********************************************
        [doubleQ,residuals] = postprocesssos(BlockedA,doubleb,doubleQ,nullThese,options);
        for constraint = 1:length(p)
            Qtemp = [];
            for i = 1:length(BlockedQ{constraint})
                Qtemp = blkdiag(Qtemp,doubleQ{constraint}{i});
            end
            Q{constraint} = Qtemp;
        end
    end
end

% *********************************************
%% EXTRACT DECOMPOSITION
% *********************************************
switch sol.problem
    case {0,1,2,3,4,5} % Well, it didn't f**k up completely at-least

        % *********************************************
        %% GENERATE MONOMIALS IN SOS-DECOMPOSITION
        % *********************************************
        for constraint = 1:length(p)

            if constraint > 1 & isequal(BlockedN{constraint},BlockedN{constraint-1}) & isequal(Blockedx{constraint},Blockedx{constraint-1}) & isequal(Blockedvarchange{constraint},Blockedvarchange{constraint-1}) & isequal(sizep(constraint),sizep(constraint-1))
                monoms{constraint} = monoms{constraint-1};
            else
                monoms{constraint} = [];
                totalN{constraint} = [];
                N = BlockedN{constraint};
                x = Blockedx{constraint};
                for i = 1:length(N)
                    % Original variable
                    for j = 1:size(N{i},1)
                        N{i}(j,:)=N{i}(j,:).*Blockedvarchange{constraint};
                    end
                    if isempty(N{i})
                        monoms{constraint} = [monoms{constraint};[]];
                    else
                        mi = kron(eye(sizep(constraint)),recovermonoms(N{i},x));
                        monoms{constraint} = [monoms{constraint};mi];
                    end
                end
                if isempty(monoms{constraint})
                    monoms{constraint}=1;
                end
            end

            % For small negative eigenvalues
            % this is a good quick'n'dirty approximation
            % Improve...use shifted eigenvalues and chol or what ever...
            if ~any(any(isnan(Q{constraint})))
                if isempty(Q{constraint})
                    Q{constraint}=0;
                    h{constraint}=0;
                else
                    usedVariables = find(any(Q{constraint},2));
                    if length(usedVariables)<length(Q{constraint})
                        Qpart = Q{constraint}(usedVariables,usedVariables);
                        [U,S,V]=svd(Qpart);
                        R = sqrt(S)*V';
                        h0 = R*monoms{constraint}(usedVariables);
                        if isa(h0,'sdpvar')
                            h{constraint} = clean(R*monoms{constraint}(usedVariables),options.sos.clean);
                            h{constraint} = h{constraint}(find(h{constraint}));
                        else
                            h{constraint} = h0;
                        end
                    else
                        [U,S,V]=svd(mid(Q{constraint}));
                        R = sqrt(S)*V';
                        h0 = R*monoms{constraint};

                        if isa(h0,'sdpvar')
                            h{constraint} = clean(R*monoms{constraint},options.sos.clean);
                            h{constraint} = h{constraint}(find(sum(h{constraint},2)),:);
                        else
                            h{constraint} = h0;
                        end
                    end
                end
                if isempty(ParametricVariables)
                    ParametricVariables = [];
                end
                setsos(p{constraint},h{constraint},ParametricVariables,Q{constraint},monoms{constraint});
            else
                if options.verbose>0;
                    if UncertainData
                        disp(' ');
                        disp('-> Only partial decomposition is returned (since you have uncertain data).');
                        disp('');
                    else
                        disp(' ');
                        disp('-> FAILURE : SOS decomposition not available.');
                        disp('-> The reason is probably that you are using a solver that does not deliver a dual (LMILAB)');
                        disp('-> Use sdsettings(''sos.model'',2) to circumvent this, or use another solver (SDPT3, SEDUMI,...)');
                        disp('');
                        disp('-> An alternative reason is that YALMIP detected infeasibility during the compilation phase.');
                    end
                end
            end
        end

        m = monoms;

    otherwise
        Q = [];
        m = [];
end

% Don't need these outside
yalmip('cleardual')

% Done with YALMIP, this is the time it took, minus solver
if ~isfield(sol,'solvertime')
    sol.solvertime = 0;
end

sol.yalmiptime = etime(clock,yalmip_time)-sol.solvertime;





function [exponent_m2,N_unique] = expandmatrix(exponent_m2,N_unique,n);

