function [F,failure,cause,ALREADY_MODELLED] = expandmodel(F,h,options,w)

% FIX : Current code experimental, complex, conservative, has issues with
% nonlinearities and is slow...
%
% If it wasn't for polynomials and sigmonials, it would be trivial, but the
% code is extremely messy since we want to have this functionality too

global LUbounds
global ALREADY_MODELLED
global MARKER_VARIABLES
global DUDE_ITS_A_GP
global ALREADY_MODELLED_INDEX
global REMOVE_THESE_IN_THE_END
global OPERATOR_IN_POLYNOM
global CONSTRAINTCUTSTATE

% All extended variables in the problem. It is expensive to extract this
% one so we will keep it and pass it along in the recursion
extendedvariables = yalmip('extvariables');

integers = yalmip('intvariables');
binaries = yalmip('binvariables');
if ~isempty(F)
    declareI = find(is(F,'integer'));
    declareB = find(is(F,'binary'));
    if ~isempty(declareI)
        integers = [integers depends(F(declareI))];
    end
    if ~isempty(declareB)
        binaries = [binaries depends(F(declareB))];
    end
end
yalmip('settempintvariables',integers);
yalmip('settempbinvariables',binaries);

% We keep track of all auxilliary variables introduced by YALMIP
nInitial = yalmip('nvars');

boundsAlreadySet = 0;
% Meta constraints are expanded first (iff, implies, alldifferent etc)
if ~isempty(F)
    meta = find(is(F,'meta'));
    if ~isempty(meta)
        LUbounds=setupBounds(F,options,extendedvariables);        
        boundsAlreadySet = 1;
        F = expandmeta(F);
    end
end

if isa(F,'constraint')
    F = lmi(F);
end

if nargin < 3
    options = sdpsettings;
end

if isempty(options)
    options = sdpsettings;
end

if isfield(options,'reusemodel')
    
else
    
    ALREADY_MODELLED = [];
    
    % Temporary hack to deal with a bug in CPLEX. For the implies operator (and
    % some more) YALMIP creates a dummy variable x with (x==1). Cplex fails
    % to solve problem with these stupid variables kept, hence we need to
    % remove these variables and constraints...
    MARKER_VARIABLES = [];
    
    % Temporary hack to deal with geometric programs. GPs are messy here,
    % becasue we can by mistake claim nonconvexity, since we may have no
    % sigmonial terms but indefinite quadratic term, but the whole problem is
    % meant to be solved using a GP solver. YES, globals suck, but this is
    % only temporary...hrm.
    DUDE_ITS_A_GP = 0;
    
    % Keep track of expressions that already have been modelled. Note that if a
    % graph-model already has been constructed but we now require a milp, for
    % numerical reasons, we should remove the old graph descriptions (important
    % for MPT models in particular)
    % FIX: Pre-parse the whole problem etc (solves the issues with GP also)
    ALREADY_MODELLED = {};
    ALREADY_MODELLED_INDEX = [];
    REMOVE_THESE_IN_THE_END = [];
    
    % Nonlinear operator variables are not allowed to be used in polynomial
    % expressions, except if they are exactly modelled, i.e. modelled using
    % MILP models. We will expand the model and collect variables that are in
    % polynomials, and check in the end if they are exaclty modelled
    OPERATOR_IN_POLYNOM = [];
end
% Assume success
failure = 0;
cause = '';

% Early bail out
if isempty(extendedvariables)
    return
end

if nargin < 4
    w = [];
end

if isempty(F) & isempty(h)
    return
end

% Check if it already has ben expanded in robustify or similar
F_alreadyexpanded = [];
if ~isempty(F)
    F_alreadyexpanded = [];
    already_expanded = expanded(F);
    if all(already_expanded)
        if isempty(setdiff(getvariables(h),expanded(h)))
            return
        end
    elseif any(already_expanded)
        F_alreadyexpanded = F(find(already_expanded));
        F = F(find(~already_expanded));
    end
end

if ~isempty(F)
    % Extract all simple bounds from the model, and update the internal bounds
    % in YALMIP. This is done in order to get tighter big-M models
    
    if boundsAlreadySet == 0;
        LUbounds = setupBounds(F,options,extendedvariables);
    end
    
    % Expand equalities first, since these might generate nonconvex models,
    % thus making it unnecessaryu to generate epigraphs etc
    equalities = is(F,'equality');
    if any(equalities)
        F = [F(find(equalities));F(find(~equalities))];
    end
end

% All variable indicies used in the problem
v1 = getvariables(F);
v2 = depends(F);
v3 = getvariables(h);
v4 = depends(h);

% HACK: Performance for LARGE-scale dualizations
if isequal(v3,v4) & isequal(v1,v2)
    variables = uniquestripped([v1 v3]);
else
    variables = uniquestripped([v1 v2 v3 v4]);
end

% Index to variables modeling operators
extended = find(ismembcYALMIP(variables,extendedvariables));

if nargin < 3
    options = sdpsettings;
end

% This is a tweak to allow epxansion of bilinear terms in robust problems,
% is expression such as abs(x*w) < 1 for all -1 < w < 1
% This field is set to 1 in robustify and tells YALMIP to skip checking for
% polynomial nonconvexity in the propagation
if ~isfield(options,'expandbilinear')
    options.expandbilinear = 0;
end

% Monomial information. Expensive to retrieve, so we pass this along
[monomtable,variabletype] = yalmip('monomtable');

% Is this trivially a GP, or meant to be at least?
if strcmpi(options.solver,'gpposy') | strcmpi(options.solver,'fmincon-geometric') | strcmpi(options.solver,'mosek-geometric')
    DUDE_ITS_A_GP = 1;
else
    if ~isequal(options.solver,'fmincon') & ~isequal(options.solver,'snopt') & ~isequal(options.solver,'ipopt') & ~isequal(options.solver,'') &  ~isequal(options.solver,'mosek')
        % User has specified some other solver, which does not
        % support GPs, hence it cannot be intended to be a GP
        DUDE_ITS_A_GP = 0;
    else
        % Check to see if there are any sigmonial terms on top-level
        DUDE_ITS_A_GP = ~isempty(find(variabletype(variables) == 4));
    end
end

% Constraints generated during recursive process to model operators
F_expand = ([]);

if isempty(F)
    F = ([]);
end

% First, check the objective
variables = uniquestripped([depends(h) getvariables(h)]);
monomtable = monomtable(:,extendedvariables);

% However, some of the variables are already expanded (expand can be called
% sequentially from solvemp and solverobust)
variables = setdiff1D(variables,expanded(h));

% Determine if we should aim for MILP/EXACT model directly
if options.allowmilp == 2
    method = 'exact';
else
    method = 'graph';
end

if DUDE_ITS_A_GP == 1
    options.allowmilp = 0;
    method = 'graph';
end

% Test for very common special case with only norm expression
ExtendedMap = yalmip('extendedmap');
fail = 0;
if  0%length(ExtendedMap) > 0 &&  all(strcmp('norm',{ExtendedMap.fcn}))
    for i = 1:length(ExtendedMap)
        if ~isequal(ExtendedMap(i).arg{2},2)
            fail = 1;
            break;
        end
        if ~isreal(ExtendedMap(i).arg{1})
            fail = 1;
            break;
        end
        if any(ismembcYALMIP(getvariables(ExtendedMap(i).arg{1}),extendedvariables))
             fail = 1;
            break;
        end
    end
    for i = 1:length(ExtendedMap)
        F_expand = [F_expand, cone(ExtendedMap(i).arg{1},ExtendedMap(i).var)];
    end
    F = F + lifted(F_expand,1);
    return
end
% *************************************************************************
% OK, looks good. Apply recursive expansion on the objective
% *************************************************************************
index_in_extended = find(ismembcYALMIP(variables,extendedvariables));
allExtStructs = yalmip('extstruct');
if ~isempty(index_in_extended)
    [F_expand,failure,cause] = expand(index_in_extended,variables,h,F_expand,extendedvariables,monomtable,variabletype,'objective',0,options,method,[],allExtStructs,w);
end

% *************************************************************************
% Continue with constraints
% *************************************************************************
constraint = 1;
all_extstruct = yalmip('extstruct');
while constraint <=length(F) & ~failure
    
    if ~already_expanded(constraint)       
        Fconstraint = F(constraint);
        variables = uniquestripped([depends(Fconstraint) getvariables(Fconstraint)]);
        
        % If constraint is a cut, all generated constraints must be marked
        % as cuts too
        CONSTRAINTCUTSTATE =  getcutflag(Fconstraint);
        [ix,jx,kx] = find(monomtable(variables,:));
        if ~isempty(jx) % Bug in 6.1
            if any(kx>1)
                OPERATOR_IN_POLYNOM = [OPERATOR_IN_POLYNOM extendedvariables(jx(find(kx>1)))];
            end
        end
        
        index_in_extended = find(ismembcYALMIP(variables,extendedvariables));
        if ~isempty(index_in_extended)
            if is(Fconstraint,'equality')
                if options.allowmilp | options.allownonconvex
                    [F_expand,failure,cause] = expand(index_in_extended,variables,-sdpvar(Fconstraint),F_expand,extendedvariables,monomtable,variabletype,['constraint #' num2str(constraint)],0,options,'exact',[],allExtStructs,w);
                else
                    failure = 1;
                    cause = ['integer model required for equality in constraint #' num2str(constraint)];
                end
            else
                [F_expand,failure,cause] = expand(index_in_extended,variables,-sdpvar(Fconstraint),F_expand,extendedvariables,monomtable,variabletype,['constraint #' num2str(constraint)],0,options,method,[],allExtStructs,w);
            end
        end
    end
    constraint = constraint+1;
end

CONSTRAINTCUTSTATE = 0;

% *************************************************************************
% Temporary hack to fix the implies operator (cplex has some problem on
% these trivial models where a variable only is used in x==1
% FIX: Automatically support this type of nonlinear operators
% *************************************************************************
if ~isempty(MARKER_VARIABLES)
    MARKER_VARIABLES = sort(MARKER_VARIABLES);
    equalities = find(is(F,'equality'));
    equalities = equalities(:)';
    remove = [];
    for j = equalities
        v = getvariables(F(j));
        if length(v)==1
            if ismembcYALMIP(v,MARKER_VARIABLES)
                remove = [remove j];
            end
        end
    end
    if ~isempty(remove)
        F(remove) = [];
    end
end

nNow = yalmip('nvars');
if nNow > nInitial
    % YALMIP has introduced auxilliary variables
    % We mark these as auxilliary
    yalmip('addauxvariables',nInitial+1:nNow);
end


F_expand = lifted(F_expand,1);
% *************************************************************************
% We are done. We might have generated some stuff more than once, but
% luckily we keep track of these mistakes and remove them in the end (this
% happens if we have constraints like (max(x)<1) + (max(x)>0) where
% the first constraint would genrate a graph-model but the second set
% creates a integer model.
% *************************************************************************
if ~failure
    F = F + F_expand;
    if length(REMOVE_THESE_IN_THE_END) > 0
        F = F(find(~ismember(getlmiid(F),REMOVE_THESE_IN_THE_END)));
    end
end

% *************************************************************************
% Normally, operators are not allowed in polynomial expressions. We do
% however allow this if the variable has been modelled with an exact MILP
% model.
% *************************************************************************
if ~failure
    dummy = unique(OPERATOR_IN_POLYNOM);
    if ~isempty(dummy)
        for i = 1:length(dummy)
            aux(i,1) = find(ALREADY_MODELLED_INDEX == dummy(i));
        end
        %    Final_model = {ALREADY_MODELLED{unique(OPERATOR_IN_POLYNOM)}};
        Final_model = {ALREADY_MODELLED{aux}};
        for i = 1:length(Final_model)
            if ~(strcmp(Final_model{i}.method,'integer') | strcmp(Final_model{i}.method,'exact') | options.allownonconvex)
                failure = 1;
                cause = 'Nonlinear operator in polynomial expression.';
                return
            end
        end
    end
end

% Append the previously appended
F = F + F_alreadyexpanded;
% Declare this model as expanded
F = expanded(F,1);

function [F_expand,failure,cause] = expand(index_in_extended,variables,expression,F_expand,extendedvariables,monomtable,variabletype,where,level,options,method,extstruct,allExtStruct,w)
global DUDE_ITS_A_GP ALREADY_MODELLED ALREADY_MODELLED_INDEX REMOVE_THESE_IN_THE_END OPERATOR_IN_POLYNOM

% *************************************************************************
% Go through all parts of expression to check for convexity/concavity
% First, a small gateway function before calling the recursive stuff
% *************************************************************************
if ~DUDE_ITS_A_GP
    [ix,jx,kx] = find(monomtable(variables,:));
    if ~isempty(jx) % Bug in 6.1
        if any(kx>1)
            OPERATOR_IN_POLYNOM = [OPERATOR_IN_POLYNOM extendedvariables(jx(find(kx>1)))];
        end
    end
end

failure = 0;
j = 1;

try
    % Optimized for 6.5 and higher (location in ismember). If user has
    % older version, it will crash and go to code below
    expression_basis = getbase(expression);
    expression_vars  = getvariables(expression);
    [yesno,location] = ismember(variables(index_in_extended),expression_vars);
    ztemp = recover(variables(index_in_extended));
    while j<=length(index_in_extended) & ~failure
        % i = index_in_extended(j);
        % zi = recover(variables(i));
        zi = ztemp(j);%recover(variables(i));
        basis = expression_basis(:,1 + location(j));
        if all(basis == 0) % The nonlinear term is inside a monomial
            if options.allownonconvex
                [F_expand,failure,cause] = expandrecursive(zi,F_expand,extendedvariables,monomtable,variabletype,where,level+1,options,'exact',[],'exact',allExtStruct,w);
            else
                failure = 1;
                cause = 'Possible nonconvexity due to operator in monomial';
            end
            %[F_expand,failure,cause] = expandrecursive(zi,F_expand,extendedvariables,monomtable,variabletype,where,level+1,options,method,[],'exact',allExtStruct,w);
        elseif all(basis >= 0)
            [F_expand,failure,cause] = expandrecursive(zi,F_expand,extendedvariables,monomtable,variabletype,where,level+1,options,method,[],'convex',allExtStruct,w);
        else
            [F_expand,failure,cause] = expandrecursive(zi,F_expand,extendedvariables,monomtable,variabletype,where,level+1,options,method,[],'concave',allExtStruct,w);
        end
        j=j+1;
    end
catch
    while j<=length(index_in_extended) & ~failure
        i = index_in_extended(j);
        basis = getbasematrix(expression,variables(i));
        if all(basis >= 0)
            [F_expand,failure,cause] = expandrecursive(recover(variables(i)),F_expand,extendedvariables,monomtable,variabletype,where,level+1,options,method,[],'convex',allExtStruct,w);
        else
            [F_expand,failure,cause] = expandrecursive(recover(variables(i)),F_expand,extendedvariables,monomtable,variabletype,where,level+1,options,method,[],'concave',allExtStruct,w);
        end
        j=j+1;
    end
end
