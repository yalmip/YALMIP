function [F,obj,m,everything] = compilesos(F,obj,options,params,candidateMonomials)
%COMPILESOS Derive sum-of-squares model without solving
%
%    [F,obj,m] = compilesos(F,h,options,params,monomials) compiles the SOS
%    problem (i.e., derives the SDP model) without actually solving it
%
%    Inputs
%     F         : The model involving SOS constraints
%     h         : Objective function (function of params) [optional]
%     options   : SDPSETTINGS structure [optional]
%     params    : Parametric variables in model [optional]
%     monomials : Prespecified monomials to be used [optional]
%
%    Outputs
%     F         : Constraints defining the problem
%     h         : Objective function
%     m         : Monomials used in the decomposition
%
% NOTE: If you use compilesos together with optimizer to solve many sos
% problems repeatedly, you must set sos.model option to 2. This is done
% automatically if you define a sos problem directly through optimizer,
% thus bypassing compilesos.  
% 
% See also OPTIMIZE, SOS, OPTIMIZER

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

if isempty(options)
    options = sdpsettings;
end

if ~isempty(options)
    if options.sos.numblkdg
        error('Does not make sense to ask for numerical block diagonalization when compiling model');
        return
    end
end

% Lazy syntax (not official...)
if nargin==1 & isa(F,'sdpvar')
    F = (sos(F));
end

% Default return structure
everything.p = [];
everything.sizep = [];
everything.normp = [];
everything.BlockedA = [];
everything.Blockedb = [];
everything.BlockedN = [];
everything.Blockedx = [];
everything.Blockedvarchange = [];
everything.BlockedQ = [];
everything.ranks = [];
everything.ParametricVariables = [];
everything.UncertainData = [];
everything.sol.problem = 0;

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
    F_parametric = ([]);
end

% *************************************************************************
%% Expand the parametric constraints
% *************************************************************************
if ~isempty(yalmip('extvariables'))
    [F_parametric,failure] = expandmodel(F_parametric,obj,options);
    F_parametric = expanded(F_parametric,1);
    obj = expanded(obj,1);    
    if failure
        error('Could not expand the model');
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
ParametricVariables = uniquestripped([depends(obj) depends(F_parametric) depends(params)]);

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
F_sos = ([]);

% *************************************************************************
%% FIGURE OUT ALL USED PARAMETRIC VARIABLES
% *************************************************************************
AllVariables =  uniquestripped([depends(obj) depends(F_original) depends(F_parametric)]);
ParametricVariables = intersect(ParametricVariables,AllVariables);
MonomVariables = setdiff(AllVariables,ParametricVariables);
params = recover(ParametricVariables);
if isempty(MonomVariables)
    error('No independent variables? Perhaps you added a constraint (p(x)) when you meant (sos(p(x))). It could also be that you added a constraint directly in the independents, such as p(x)>=0 or similarily.');
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
options = selectSOSmodel(F,options,NonLinearParameterization,noRANK,IntegerData,UncertainData);

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
    tempops.usex0 = 0;
    [aux1,aux2,aux3,LPmodel] = export((temp>=0),temp,tempops);   
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
        % This is the case (sos(t)) where t is a parametric (matrix) variable
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
        %  error('You have constraints of the type (sos(f(parametric_variables))). Please use (f(parametric_variables) > 0) instead')
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
            sol.yalmiptime = 0;
            sol.solvertime = 0;
            sol.info = yalmiperror(1,'YALMIP');
            sol.problem = 2;
            everything.sol = sol;
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

    if options.sos.reuse & constraint > 1 && isequal(previous_exponent_p_monoms,exponent_p_monoms)
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
BlockedQ = [];
switch options.sos.model
    case 1
        % Kernel model
        [F,obj,BlockedQ,Primal_matrices,Free_variables] = create_kernelmodel(BlockedA,Blockedb,F_parametric,parobj,options,[]);
    case {2,4,5,6}
        % 2=Image model, 4=reduced nonlinear, 5=dd,6=sd
        [F,obj,BlockedQ,sol] = create_imagemodel(BlockedA,Blockedb,F_parametric,parobj,options);

    case 3
        % Un-official model to solve bilinearly parameterized SOS using SDPLR
        [F,obj,options] = create_lrmodel(BlockedA,Blockedb,F_parametric,parobj,options,ParametricVariables);

    otherwise
end

for constraint = 1:length(p)
    if constraint > 1 && isequal(BlockedN{constraint},BlockedN{constraint-1}) && isequal(Blockedx{constraint},Blockedx{constraint-1}) && isequal(Blockedvarchange{constraint},Blockedvarchange{constraint-1}) && isequal(sizep(constraint),sizep(constraint-1))
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

end
m = monoms;

everything.p = p;
everything.sizep = sizep;
everything.normp = normp;
everything.BlockedA = BlockedA;
everything.Blockedb = Blockedb;
everything.BlockedN = BlockedN;
everything.Blockedx = Blockedx;
everything.Blockedvarchange = Blockedvarchange;
everything.BlockedQ = BlockedQ;
everything.ranks = ranks;
everything.ParametricVariables = ParametricVariables;
everything.UncertainData = UncertainData;
everything.sol = sol;