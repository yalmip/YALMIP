function [F,obj,m] = solvesos(F,obj,options,params,candidateMonomials)
%COMPILESOS Sum of squares decomposition
%
%    [F,obj,m] = compilesos(F,h,options,params,monomials) derives the SOS
%    problem without actually solving it
% 
% See also SOLVESOS

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
ParametricVariables = uniquestripped([depends(obj) depends(F_parametric) depends(params)]);

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
        noCOMPLICATING = ~any(ismember([7 8 9 10 12 13 14 15],constraint_classes));
        if noCOMPLICATING & ~NonLinearParameterization & noRANK & ~IntegerData
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
        %% PRUNE W.R.T USER DEFINED
        %*********************************************
        %         if ~isempty(exponent_c)
        %             hash = randn(size(exponent_m{1},2),1);
        %             for i = 1:length(exponent_m)
        %                 hash = randn(size(exponent_m{i},2),1);
        %                 keep = find(ismember(exponent_m{i}*hash,exponent_c{1}*hash));
        %                 if length(keep)>0
        %                     exponent_m{i} = exponent_m{i}(keep,:);
        %                 else
        %                     exponent_m{i}=[];
        %                 end
        %             end
        %         end

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
    case {2,4}
        % Image model
        [F,obj,BlockedQ,sol] = create_imagemodel(BlockedA,Blockedb,F_parametric,parobj,options);

    case 3
        % Un-official model to solve bilinearly parameterized SOS using SDPLR
        [F,obj,options] = create_lrmodel(BlockedA,Blockedb,F_parametric,parobj,options,ParametricVariables);

    otherwise
end

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

end
m = monoms;





function p_base_parametric = parameterizedbase(p,z, params,ParametricIndicies,exponent_p,p_base)

% Check for linear parameterization
parametric_basis = exponent_p(:,ParametricIndicies);
if all(sum(parametric_basis,2)==0)
    p_base_parametric = full(p_base(:));
    return
end
if all(sum(parametric_basis,2)<=1)
    p_base_parametric = full(p_base(:));
    n = length(p_base_parametric);


    if 1
        [ii,vars] = find(parametric_basis);
        ii = ii(:)';
        vars = vars(:)';
    else
        ii = [];
        vars = [];
        js = sum(parametric_basis,1);
        indicies = find(js);
        for i = indicies
            if js(i)
                j = find(parametric_basis(:,i));
                ii = [ii j(:)'];
                vars = [vars repmat(i,1,js(i))];
            end
        end
    end

    k = setdiff1D(1:n,ii);
    if isempty(k)
        p_base_parametric = p_base_parametric.*sparse(ii,repmat(1,1,n),params(vars));
    else
        pp = params(vars); % Must do this, bug in ML 6.1 (x=sparse(1);x([1 1]) gives different result in 6.1 and 7.0!)
        p_base_parametric = p_base_parametric.*sparse([ii k(:)'],repmat(1,1,n),[pp(:)' ones(1,1,length(k))]);
    end
else
    % Bummer, nonlinear parameterization sucks...
    for i = 1:length(p_base)
        j = find(exponent_p(i,ParametricIndicies));
        if ~isempty(j)
            temp = p_base(i);

            for k = 1:length(j)
                if exponent_p(i,ParametricIndicies(j(k)))==1
                    temp = temp*params(j(k));
                else
                    temp = temp*params(j(k))^exponent_p(i,ParametricIndicies(j(k)));
                end
            end
            xx{i} = temp;
        else
            xx{i} = p_base(i);
        end
    end
    p_base_parametric = stackcell(sdpvar(1,1),xx)';
end



function [A,b] = generate_kernel_representation_data(N,N_unique,exponent_m2,exponent_p,p,options,p_base_parametric,ParametricIndicies,MonomIndicies,FirstRun)

persistent saveData

exponent_p_parametric = exponent_p(:,ParametricIndicies);
exponent_p_monoms = exponent_p(:,MonomIndicies);
pcoeffs = getbase(p);
if any(exponent_p_monoms(1,:))
    pcoeffs=pcoeffs(:,2:end); % No constant term in p
end
b = [];

parametric = full((~isempty(ParametricIndicies) & any(any(exponent_p_parametric))));

% For problems with a lot of similar cones, this saves some time
reuse = 0;
if ~isempty(saveData) & isequal(saveData.N,N) & ~FirstRun
    n = saveData.n;
    ind = saveData.ind;
    if  isequal(saveData.N_unique,N_unique) & isequal(saveData.exponent_m2,exponent_m2)% & isequal(saveData.epm,exponent_p_monoms)
        reuse = 1;
    end
else
    % Congruence partition sizes
    for k = 1:size(N,1)
        n(k) = size(N{k},1);
    end
    % Save old SOS definition
    saveData.N = N;
    saveData.n = n;
    saveData.N_unique = N_unique;
    saveData.exponent_m2 = exponent_m2;
    saveData.N_unique = N_unique;
end

if reuse & options.sos.reuse
    % Get old stuff
    if size(exponent_m2{1},2)==2 % Stupid set(sos(parametric)) case
        ind = spalloc(1,1,0);
        ind(1)=1;
        allj = 1:size(exponent_p_monoms,1);
        used_in_p = ones(size(exponent_p_monoms,1),1);
    else
        allj = [];
        used_in_p = zeros(size(exponent_p_monoms,1),1);
        hash = randn(size(exponent_p_monoms,2),1);
        exponent_p_monoms_hash = exponent_p_monoms*hash;
        for i = 1:size(N_unique,1)
            monom = sparse(N_unique(i,3:end));
            j = find(exponent_p_monoms_hash == (monom*hash));
          
            if isempty(j)
                b = [b 0];
                allj(end+1,1) = 0;
            else
                used_in_p(j) = 1;
                allj(end+1,1:length(j)) = j(:)';
            end
        end
        ind = saveData.ind;
    end
else
    allj = [];
    used_in_p = zeros(size(exponent_p_monoms,1),1);
    if size(exponent_m2{1},2)==2 % Stupid set(sos(parametric)) case
        ind = spalloc(1,1,0);
        ind(1)=1;
        allj = 1:size(exponent_p_monoms,1);
        used_in_p = ones(size(exponent_p_monoms,1),1);
    else
        % To speed up some searching, we random-hash data
        hash = randn(size(exponent_p_monoms,2),1);
        for k = 1:length(exponent_m2)
            if isempty(exponent_m2{k})
                exp_hash{k}=[];
            else
                exp_hash{k} = sparse((exponent_m2{k}(:,3:end)))*hash; % SPARSE NEEDED DUE TO STRANGE NUMERICS IN MATLAB ON 0s (the stuff will differ on last bit in hex format)
            end
        end

        p_hash = exponent_p_monoms*hash;
        ind = spalloc(size(N_unique,1),sum(n.^2),0);       

        for i = 1:size(N_unique,1)
            monom = N_unique(i,3:end);
            monom_hash = sparse(monom)*hash;
            LHS = 0;
            start = 0;
            for k = 1:size(N,1)
                j = find(exp_hash{k} == monom_hash);
                if ~isempty(j)
                    pos=exponent_m2{k}(j,1:2);
                    nss = pos(:,1);
                    mss = pos(:,2);
                    indicies = nss+(mss-1)*n(k);
                    ind(i,indicies+start) = ind(i,indicies+start) + 1;                         
                end
                start = start + (n(k))^2;
                %                start = start + (matrixSOSsize*n(k))^2;
            end

            j = find(p_hash == monom_hash);

            if isempty(j)
                allj(end+1,1) = 0;
            else
                used_in_p(j) = 1;
                allj(end+1,1:length(j)) = j(:)';
            end
        end
    end
end
saveData.ind = ind;

% Some parametric terms in p(x,t) do not appear in v'Qv
% So these have to be added 0*Q = b
not_dealt_with  = find(used_in_p==0);
while ~isempty(not_dealt_with)
    j = findrows(exponent_p_monoms,exponent_p_monoms(not_dealt_with(1),:));
    allj(end+1,1:length(j)) = j(:)';
    used_in_p(j) = 1;
    not_dealt_with  = find(used_in_p==0);
    ind(end+1,1)=0;
end

matrixSOSsize = length(p);
if parametric
    % Inconsistent behaviour in MATLAB
    if size(allj,1)==1
        uu = [0;p_base_parametric];
        b = sum(uu(allj+1))';
    else
        b = [];
        for i = 1:matrixSOSsize
            for j = i:matrixSOSsize
                if i~=j
                    uu = [0;2*p_base_parametric(:,(i-1)*matrixSOSsize+j)];
                else
                    uu = [0;p_base_parametric(:,(i-1)*matrixSOSsize+j)];
                end
                b = [b sum(uu(allj+1),2)'];
            end
        end
    end
else
    if matrixSOSsize == 1
        uu = [zeros(size(pcoeffs,1),1) pcoeffs]';
        b = sum(uu(allj+1,:),2)';
    else
        b = [];
        for i = 1:matrixSOSsize
            for j = i:matrixSOSsize
                if i~=j
                    uu = [0;2*pcoeffs((i-1)*matrixSOSsize+j,:)'];
                else
                    uu = [0;pcoeffs((i-1)*matrixSOSsize+j,:)'];
                end
                b = [b;sum(uu(allj+1,:),2)'];
            end
        end
    end
    % uu = [0;pcoeffs(:)];
    % b = sum(uu(allj+1),2)';
end

b = b';
dualbase = ind;

j = 1;
A = cell(size(N,1),1);
for k = 1:size(N,1)
    if matrixSOSsize==1
        A{k} = dualbase(:,j:j+n(k)^2-1);
    else
        % Quick fix for matrix SOS case, should be optimized
        A{k} = inflate(dualbase(:,j:j+n(k)^2-1),matrixSOSsize,n(k));
    end
    j = j + n(k)^2;
end
b = b(:);



function newAi = inflate(Ai,matrixSOSsize,n);
% Quick fix for matrix SOS case, should be optimized
newAi = [];
for i = 1:matrixSOSsize
    for r = i:matrixSOSsize
        for m = 1:size(Ai,1)
            ai = reshape(Ai(m,:),n,n);
            V = zeros(matrixSOSsize,matrixSOSsize);
            V(i,r)=1;
            V(r,i)=1;
            ai = kron(V,ai);
            newAi = [newAi;ai(:)'];
        end
    end
end
        

function [Z,Q1,R] = sparsenull(A)

[Q,R] = qr(A');
n = max(find(sum(abs(R),2)));
Q1 = Q(:,1:n);
R = R(1:n,:);
Z = Q(:,n+1:end); % New basis


function [F,obj,BlockedQ,Primal_matrices,Free_variables] = create_kernelmodel(BlockedA,Blockedb,F_parametric,parobj,options,sparsityPattern);

% To get the primal kernel representation, we simply use
% the built-in dualization module!
% First, write the problem in primal kernel format
traceobj = 0;
dotraceobj = options.sos.traceobj;
F = F_parametric;
for i = 1:length(Blockedb)


    sizematrixSOS = sqrt(size(Blockedb{i},2));
    for k = 1:sizematrixSOS
        for r = k:sizematrixSOS
            res{(k-1)*sizematrixSOS+r} = 0;
        end
    end

    for j = 1:length(BlockedA{i})
        n = sqrt(size(BlockedA{i}{j},2));
        BlockedQ{i}{j} = sdpvar(n*sizematrixSOS,n*sizematrixSOS);
        F = F + set(BlockedQ{i}{j});
        if sizematrixSOS>0
            % Matrix valued sum of sqaures
            % Loop over all elements
            starttop = 1;
            for k = 1:sizematrixSOS
                startleft = 1;
                for r = 1:sizematrixSOS
                    if k<=r
                        Qkr = BlockedQ{i}{j}(starttop:starttop+n-1,startleft:startleft+n-1);
                        res{(k-1)*sizematrixSOS+r} = res{(k-1)*sizematrixSOS+r} + BlockedA{i}{j}*reshape(Qkr,n^2,1);
                    end
                    startleft = startleft + n;
                end
                starttop = starttop + n;
            end
        else
            % Standard case
            res{1} = res{1} + BlockedA{i}{j}*reshape(BlockedQ{i}{j},n^2,1);
        end
        if dotraceobj
            traceobj = traceobj + trace(BlockedQ{i}{j});
        end
    end
    for k = 1:sizematrixSOS
        for r = k:sizematrixSOS
            F = F + set(res{(k-1)*sizematrixSOS+r} == Blockedb{i}(:,(k-1)*sizematrixSOS+r));
        end
    end
end

% % Constrain elements according to a desired sparsity
if ~isempty(sparsityPattern)
    res = [];
    for i = 1:length(BlockedQ)
        for j = 1:length(BlockedQ{i})
            if ~isempty(sparsityPattern{i}{j})
                H = spalloc(length(BlockedQ{i}{j}),length(BlockedQ{i}{j}),length(sparsityPattern{i}{j}));
                H(sparsityPattern{i}{j}) = 1;
                k = find(triu(H));
                res = [res;BlockedQ{i}{j}(k)];
            end
        end
    end
    F = F + set(0 == res);
end

% And get the primal model of this
if isempty(parobj)
    if options.sos.traceobj
        [F,obj,Primal_matrices,Free_variables] = dualize(F,traceobj,1,options.sos.extlp);
    else
        [F,obj,Primal_matrices,Free_variables] = dualize(F,[],1,options.sos.extlp);
    end
else
    [F,obj,Primal_matrices,Free_variables] = dualize(F,parobj,1,options.sos.extlp);
end
% In dual mode, we maximize
obj = -obj;


function [F,obj,BlockedQ,sol] = create_imagemodel(BlockedA,Blockedb,F_parametric,parobj,options);


% *********************************
% IMAGE REPRESENTATION
% Needed for nonlinearly parameterized problems
% More efficient in some cases
% *********************************
g = [];
vecQ = [];
sol.problem = 0;
for i = 1:length(BlockedA)
    q = [];
    A = [];
    for j = 1:length(BlockedA{i})
        n = sqrt(size(BlockedA{i}{j},2));
        BlockedQ{i}{j} = sdpvar(n,n);
        q = [q;reshape(BlockedQ{i}{j},n^2,1)];
        A = [A BlockedA{i}{j}];
    end
    vecQ = [vecQ;q];
    g = [g;Blockedb{i}-A*q];
end

g_vars = getvariables(g);
q_vars = getvariables(vecQ);
x_vars = setdiff(g_vars,q_vars);

base = getbase(g);
if isempty(x_vars)
    A = base(:,1);base = base(:,2:end);
    B = (base(:,ismember(g_vars,q_vars)));
    Bnull = sparsenull(B);
    t = sdpvar(size(Bnull,2),1);
    imQ = -B\A+Bnull*t;
else
    A = base(:,1);base = base(:,2:end);
    C = base(:,ismember(g_vars,x_vars));
    B = (base(:,ismember(g_vars,q_vars)));
    [Bnull,Q1,R1] = sparsenull(B);Bnull(abs(Bnull) < 1e-12) = 0;
    t = sdpvar(size(Bnull,2),1);
    if options.sos.model == 2
        imQ = -Q1*(R1'\(A+C*recover(x_vars)))+Bnull*t;
    else
        H = testnice(B);
        imQ = -B\(A+C*recover(x_vars))+H*t;        
    end
end
notUsed = find(sum(abs(B),2)==0);
if ~isempty(notUsed)
    ff=g(notUsed);
    if isa(ff,'double')
        if nnz(ff)>0
            sol.yalmiptime = 0; % FIX
            sol.solvertime = 0;
            sol.problem = 2;
            sol.info = yalmiperror(1,'YALMIP');
            F = [];
            obj = [];
            if options.verbose > 0
                disp(' ');
                disp('-> During construction of data, I encountered a situation')
                disp('-> situation that tells me that the problem is trivially infeasible!')
                disp('-> Have you forgotten to define some parametric variables?,')
                disp('-> or perhaps you have a parametric problem where the highest')
                disp('-> power in some of the independent variables is odd for any choice')
                disp('-> of parametric variables, such as x^8+x^7+t*y^3')
                disp('-> Anyway, take a good look at your model again...');                
            end
            return
            %            error('You seem to have a strange model. Have you forgotten to define some parametric variable?');
        end
    else
       F_parametric = F_parametric + set(g(notUsed)==0);
    end
end
F_sos = set([]);
obj = 0;
for i = 1:length(BlockedQ)
    for j = 1:size(BlockedQ{i},2)
        Q_old = BlockedQ{i}{j};
        Q_old_vars = getvariables(Q_old);
        Q_old_base = getbase(Q_old);
        in_this = [];
        for k = 1:length(Q_old_vars)
            in_this = [in_this;find(Q_old_vars(k)==q_vars)];
        end
        Q_new = Q_old_base(:,1) + Q_old_base(:,2:end)*imQ(in_this);
        Q_new = reshape(Q_new,length(BlockedQ{i}{j}),length(BlockedQ{i}{j}));
        obj = obj+trace(Q_new);
        if ~isa(Q_new,'double')
            F_sos = F_sos + set(Q_new);
        elseif min(eig(Q_new))<-1e-8
            sol.yalmiptime = 0; % FIX
            sol.solvertime = 0;
            sol.problem = 2;
            sol.info = yalmiperror(1,'YALMIP');
            F = [];
            obj = [];
            error('Problem is trivially infeasible. After block-diagonalizing, I found constant negative definite blocks!');
            return
        end
        BlockedQ{i}{j} = Q_new;
    end
end

F = F_parametric + F_sos;

if isa(obj,'double') | (options.sos.traceobj == 0)
    obj = [];
end

if ~isempty(parobj)
    obj = parobj;
end


function   [F,obj,options] = create_lrmodel(BlockedA,Blockedb,F_parametric,parobj,options,ParametricVariables)
% Some special code for ther low-rank model in SDPLR
% Experimental code, not official yet
allb = [];
allA = [];
K.s = [];
for i = 1:length(Blockedb)
    allb = [allb;Blockedb{i}];
    Ai = [];
    for j = 1:size(BlockedA{i},2)
        Ai = [Ai BlockedA{i}{j}];
        K.s = [K.s sqrt(size(BlockedA{i}{j},2))];
    end
    %blkdiag bug in 7.0...
    [n1,m1] = size(allA);
    [n2,m2] = size(Ai);
    allA = [allA spalloc(n1,m2,0);spalloc(n2,m1,0) Ai];
end
options.solver = 'sdplr';
z = recover(ParametricVariables)
start = size(BlockedA,2)+1;
Mi = [];
for i = 1:length(allb)
    if isa(allb(i),'sdpvar')
        [Qi,ci,fi,xi,infoi] = quaddecomp(allb(i),z);
    else
        Qi = zeros(length(z));
        ci = zeros(length(z),1);
        fi = allb(i);
    end
    Z = -[fi ci'/2;ci/2 Qi];
    Mi = [Mi;Z(:)'];
end
K.s = [K.s length(z)+1];
zeroRow = zeros(1,size(allA,2));
allA = [allA Mi;zeroRow 1 zeros(1,K.s(end)^2-1)];
b = zeros(size(allA,1),1);b(end) = 1;
y = sdpvar(length(b),1);
CminusAy = -allA'*y;
start = 1;

% Get the cost, expressed in Z
[Qi,ci,fi,xi,infoi] = quaddecomp(parobj,z);
C = [fi ci'/2;ci/2 Qi];
F = set([]);
for i = 1:length(K.s)
    if i<length(K.s)
        F = F + set(reshape(CminusAy(start:start+K.s(i)^2-1),K.s(i),K.s(i)));
    else
        F = F + set(reshape(C(:) + CminusAy(start:start+K.s(i)^2-1),K.s(i),K.s(i)));
    end
    start = start + K.s(i)^2;
end
obj = -b'*y;

options.sdplr.forcerank = ceil(K.s/2);
options.sdplr.forcerank(end) = 1;
options.sdplr.feastol = 1e-7;



function [exponent_m2,N_unique] = expandmatrix(exponent_m2,N_unique,n);


function [sol,m,Q,residuals,everything] = solvesos_find_blocks(F,obj,options,params,candidateMonomials)

tol = options.sos.numblkdg;
if tol > 1e-2
    disp(' ');
    disp('-> Are you sure you meant to have a tolerance in numblk that big!')
    disp('-> The options numblkdiag controls the tolerance, it is not a 0/1 switch.')
    disp(' ');
end
options.sos.numblkdg = 0;
[sol,m,Q,residuals,everything] = solvesos(F,obj,options,params,candidateMonomials);

% Save old structure to find out when we have stalled
for i = 1:length(Q)
    oldlengths{i} = length(Q{i});
end

go_on = (sol.problem == 0 | sol.problem == 4);
while go_on

    for sosfun = 1:length(Q)
        Qtemp = Q{sosfun};
        keep = diag(Qtemp)>tol;
        Qtemp(:,find(~keep)) = [];
        Qtemp(find(~keep),:) = [];

        m{sosfun} = m{sosfun}(find(keep));

        Qtemp(abs(Qtemp) < tol) = 0;
        [v1,dummy1,r1,dummy3]=dmperm(Qtemp+eye(length(Qtemp)));
        lengths{sosfun} = [];
        n{sosfun} = {};
        for blocks = 1:length(r1)-1
            i1 = r1(blocks);
            i2 = r1(blocks+1)-1;
            if i2>i1
                n{sosfun}{blocks} = m{sosfun}(v1(i1:i2));
            else
                n{sosfun}{blocks} = m{sosfun}(v1(i1));
            end
            lengths{sosfun} =  [lengths{sosfun}  length(n{sosfun}{blocks})];
        end
        lengths{sosfun} = sort(lengths{sosfun});
    end
    go_on = ~isequal(lengths,oldlengths);
    oldlengths = lengths;
    if go_on
        [sol,m,Q,residuals,everything] = solvesos(F,obj,options,params,n);
        go_on = go_on & (sol.problem == 0 | sol.problem == 4);
        if sol.problem == 1
            disp('-> Feasibility was lost during the numerical block-diagonalization.')
            disp('-> The setting sos.numblkdiag is probably too big')
        end
    end
end


function sol =  solveranksos(F,obj,options,ranks,BlockedQ)

Frank = set([]);
for i = 1:length(ranks)
    if ~isinf(ranks(i))
        Frank = Frank + set(rank(BlockedQ{i}{1}) <= ranks(i));
    end
end
% rank adds the pos.def constraints again!!, so we remove them
check = ones(length(F),1);
keep  = ones(length(F),1);
for i = 1:length(BlockedQ)
    for j = 1:length(BlockedQ{i})
        Qij = BlockedQ{i}{j};
        for k = find(check)'
            if isequal(Qij,sdpvar(F(k)))
                keep(k)  = 0;
                check(k) = 0;
            end
        end
    end
end
% Let's hope LMIRANK is there
sol =  solvesdp(F(find(keep)) + Frank,[],sdpsettings(options,'solver',''));


function  [indx]=colspaces(A)
indx = [];
for i = 1:size(A,2)
    s = max(find(A(:,i)));
    indx = [indx s];
end
indx = unique(indx);

function H = testnice(A_equ)
[L,U,P] = lu(A_equ);
[L,U,P] = lu(A_equ');
r = colspaces(L');
AA = L';
H1 = AA(:,r);
H2 = AA(:,setdiff(1:size(AA,2),r));
H = P'*[-H1\H2;eye(size(H2,2))];
    
