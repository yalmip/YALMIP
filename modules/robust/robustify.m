function [F,h,failure] = robustify(F,h,ops,w)
%ROBUSTIFY  Derives robust counterpart.
%
% [Frobust,objrobust,failure] = ROBUSTIFY(F,h,options) is used to derive
% the robust counterpart of an uncertain YALMIP model.
%
%   min        h(x,w)
%   subject to
%           F(x,w) >(=) 0  for all w in W
%
% The constraints and objective have to satisfy a number of conditions for
% the robustification to be possible. Please refer to the YALMIP Wiki for
% the current assumptions.
%
% Some options for the robustification strategies can be altered via the
% solver tag 'robust' in sdpsettings
%
%  'robust.lplp'  : Controls how linear constraints with affine
%                   parameterization in an uncertainty with polytopic
%                   description is handled. Can be either 'duality' or
%                   'enumeration' 
%
%  'robust.auxred': Controls how uncertainty dependent auxiliary variables
%                   are handled
%                   Can be either 'projection' or 'enumeration' (exact),
%                   or 'none' or 'affine' (conservative)
%
%  'robust.reducedual' Controls if the system equality constraints derived
%                   when using the duality filter should be eliminated,
%                   thus reducing the number of variables, possibly
%                   destroying sparsity .
%
%  'robust.polya'  : Controls the relaxation order of polynomials. If set to
%                   NAN, the polynomials will be eliminated by forcing the
%                   coefficients to zero
%
% See also UNCERTAIN

% Author Johan Löfberg
% $Id: robustify.m,v 1.55 2010-03-10 15:19:05 joloef Exp $

failure = 0;
if nargin < 3
    ops = sdpsettings;
elseif isempty(ops)
    ops = sdpsettings;
end

if nargin < 4
    w = [];
end

if nargin>1
    if isa(h,'double')
        h = [];
    end
else
    h = [];
end

% We keep track of auxilliary generated variables
nInitial = yalmip('nvars');

% Find the scenario, extract uncertainty model and classifiy variables
[UncertainModel,Uncertainty,VariableType,ops] = decomposeUncertain(F,h,w,ops);
x = VariableType.x;
w = VariableType.w;

if isempty(x)
    error('There are no decision variables in the uncertain model.')
end

if isempty(UncertainModel.F_xw)
    error('The uncertainty does not enter the model anywhere.');
end

% Experimental code for conic-conic case
if ops.robust.coniclp.useconicconic || ((any(is(UncertainModel.F_xw,'sdp')) ||  any(is(UncertainModel.F_xw,'socp'))) && (any(is(Uncertainty.F_w,'sdp')) ||  any(is(Uncertainty.F_w,'socp'))))
    SOSModel = [];
    for i = 1:length(UncertainModel.F_xw)
        if any(ismember(depends(UncertainModel.F_xw(i)),getvariables(VariableType.w)))
            SOSModel = [SOSModel, dualtososrobustness(UncertainModel.F_xw(i),Uncertainty.F_w,VariableType.w,VariableType.x,ops.robust.conicconic.tau_degree,ops.robust.conicconic.gamma_degree,ops.robust.conicconic.Z_degree)];
        else
            % Misplaced?
            SOSModel = [SOSModel, UncertainModel.F_xw(i)];
        end
    end
    %SOSModel = expanded(SOSModel,1);
    F = [SOSModel, UncertainModel.F_x];
    h = UncertainModel.h;
    h = expanded(h,1);
    F = expanded(F,1); % This is actually done already in expandmodel
   % h = expanded(h,1); % But this one has to be done manually
  
    return
end

% FIXME: SYNC with expandmodel?
if ~isempty(UncertainModel.F_x)
    nv = yalmip('nvars');
    yalmip('setbounds',1:nv,repmat(-inf,nv,1),repmat(inf,nv,1));
    LU = getbounds(UncertainModel.F_x);
    yalmip('setbounds',1:nv,LU(:,1),LU(:,2));
end
        
% At this point, we have to decide on the algorithm we should use for
% robustifying the constraints. There are a couple of alternatives,
% depending on uncertainty and constraints
% 1. Polya:       Polynomial uncertainty dependence, simplex uncertainty,
%                 can only be applied on LP constraints
% 2. Elimination: Last resort, tries to cancel all nonlinear uncertainties
%                 by setting coefficients to zero
% 3. Explicit:    Linear uncertainty dependence, box-model uncertainty, can
%                 only be applied on LP constraints
% 4. Enumeration: Linear uncertainty dependence, polytopic uncertainty,
%                 arbitrary type of constraints (convex)
% 5. Duality:     Linear uncertainty dependence, conic uncertainty, can
%                 only be applied on LP constraints
% 6. S-procedure  Special case, quadratic dependence in elementwise, one
%                 quadratic constraint in W (obsolete)
% 7. Conic conic  Subsumes S-procedure


% Robust model
F_robust = ([]);

% We begin by checking to see if the user wants to apply Polyas theorem.
% If that is the case, search for simplex structures, and apply Polyas.
if ~isnan(ops.robust.polya) & any(strcmp(Uncertainty.uncertaintyTypes,'simplex')) & ~ops.robust.forced_enumeration
    F_polya = [];
    % Recursively apply Polya relaxation w.r.t each simplex
    for i = find(strcmp(Uncertainty.uncertaintyTypes,'simplex'))
        [UncertainModel.F_xw, F_polya] = filter_polya(UncertainModel.F_xw+F_polya,w(Uncertainty.uncertaintyGroups{i}),ops.robust.polya);
    end
    [UncertainModel.F_xw,F_robust] = pruneCertain(F_polya,F_robust,UncertainModel.F_xw,w);
end

% LP constraints with quadratic dependence and quadratic uncertainty region
% can be handled tightly using the S-procedure
if (all(strcmp(Uncertainty.uncertaintyTypes,'2-norm')) | all(strcmp(Uncertainty.uncertaintyTypes,'quadratic'))) & length(Uncertainty.uncertaintyTypes)==1 & ~ops.robust.forced_enumeration
    [UncertainModel.F_xw,F_sprocedure] = filter_sprocedure(UncertainModel.F_xw,w,Uncertainty.separatedZmodel,ops);
    F_robust = F_robust + F_sprocedure;
end

% There might still be nonlinearities left in the model. These have to be
% removed. We remove all terms with w-degree larger than 1
[UncertainModel.F_xw,F_elimination] = filter_eliminatation(UncertainModel.F_xw,w,1,ops);
F_robust = F_robust + F_elimination;

% Equality constraints cannot be part of an uncertain problem. Any
% dependence w.r.t w in equalities has to be removed
F_eq = extractConstraints(UncertainModel.F_xw,'equality');
UncertainModel.F_xw = UncertainModel.F_xw - F_eq;
[F_eq_left,F_eliminate_equality] = filter_eliminatation(F_eq,w,0,ops);
F_robust = F_robust + F_eliminate_equality + F_eq_left;

% The problem should now be linear in the uncertainty, with no uncertain
% equality constraints. Hence, now we apply explicit maximization,
% enumeration or duality-based robustification.

% We start with the norm balls
if ~ops.robust.forced_enumeration
    for i = 1:length(Uncertainty.uncertaintyTypes)
        if ismember(Uncertainty.uncertaintyTypes{i},{'1-norm','2-norm','inf-norm'})
            F_lp = extractConstraints(UncertainModel.F_xw,'elementwise');
            UncertainModel.F_xw = UncertainModel.F_xw - F_lp;
            F_flt = filter_normball(F_lp,Uncertainty.separatedZmodel{i},x,w(Uncertainty.uncertaintyGroups{i}),w,Uncertainty.uncertaintyTypes{i},ops,VariableType);
            [UncertainModel.F_xw,F_robust] = pruneCertain(F_flt,F_robust,UncertainModel.F_xw,w);
        end
    end
end

% Pick out the uncertain linear equalities and robustify using duality if
% user has opted for this or the uncertainty is conic.
conic = ~isequal(Uncertainty.Zmodel.K.s,0) | ~isequal(Uncertainty.Zmodel.K.q,0);
if (conic | isequal(ops.robust.lplp,'duality')) & ~ops.robust.forced_enumeration
    F_lp = extractConstraints(UncertainModel.F_xw,'elementwise');
    UncertainModel.F_xw = UncertainModel.F_xw - F_lp;
    nv = yalmip('nvars');
    F_filter = filter_duality(F_lp,Uncertainty.Zmodel,x,w,ops);
    F_robust = F_robust + F_filter;
    if isa(F_filter,'lmi') & ops.verbose
        newvars = nnz(getvariables(F_filter)>nv);      
        disp([' - Duality introduced ' num2str(newvars) ' variables, ' num2str(nnz(is(F_filter,'equality'))) ' equalities, ' num2str(nnz(is(F_filter,'elementwise'))) ' LP inqualities and ' num2str(nnz(is(F_filter,'sdp'))+nnz(is(F_filter,'socp'))) ' conic constraints']);       
    end
end

% Robustify remaining uncertain LP/SOCP/SDP constraints and robustify by
% enumeration.
F_conic = extractConstraints(UncertainModel.F_xw,{'sdp','socc','elementwise'});
UncertainModel.F_xw = UncertainModel.F_xw - F_conic;
[F_temp,enumerationfailed] = filter_enumeration(F_conic,Uncertainty.Zmodel,x,w,ops,Uncertainty.uncertaintyTypes,Uncertainty.separatedZmodel,VariableType);
if enumerationfailed
    % Reset to previous state
    UncertainModel.F_xw = UncertainModel.F_xw + F_conic;
else
    F_robust = F_robust + F_temp;
end

if enumerationfailed
    % Enumeration failed, probably due to lack of MPT. If problem is conic,
    % we are in trouble. If simple LP, we can resort to duality approach
    if conic
        if ops.verbose
            disp(' - Enumeration of uncertainty polytope failed. Missing Multiparametric Toolbox?')
        end
        error('Enumeration failed (lacking MPT?),  and due to conic constraints, duality cannot be used');
    else
        F_lp = extractConstraints(UncertainModel.F_xw,'elementwise');
        UncertainModel.F_xw = UncertainModel.F_xw - F_lp;
        nv = yalmip('nvars');
        F_filter = filter_duality(F_lp,Uncertainty.Zmodel,x,w,ops);
        if ops.verbose
            if isa(F_filter,'lmi')
             disp([' - Duality introduced ' num2str(yalmip('nvars')-nv') ' variables, ' num2str(nnz(is(F_filter,'equality'))) ' equalities, ' num2str(nnz(is(F_filter,'elementwise'))) ' LP inqualities and ' num2str(nnz(is(F_filter,'sdp'))+nnz(is(F_filter,'socp'))) ' conic constraints']);
            end
        end
        F_robust = F_robust + F_filter;
    end
end

% If there is anything left now, it means that we do not support it (such
% as conic uncertainty in conic constraint)
if length(UncertainModel.F_xw) > 0
    if any(~islinear(UncertainModel.F_xw))
        error('There are some uncertain constraints which cannot be robustified by YALMIP')
    else
        F_robust = F_robust + UncertainModel.F_xw;
    end
end

% Return the robustfied model
F = F_robust+UncertainModel.F_x;
h = UncertainModel.h;

% The model has been expanded, so we have to remember this (trying to
% expand an expanded model leads to nonconvexity error)
F = expanded(F,1); % This is actually done already in expandmodel
h = expanded(h,1); % But this one has to be done manually

nNow = yalmip('nvars');
if nNow > nInitial
    % YALMIP has introduced auxilliary variables
    % We mark these as auxilliary
    yalmip('addauxvariables',nInitial+1:nNow);
end

if ops.verbose
    disp('***** Derivation of robust counterpart done ***********************');
end


function [F_xw,F_robust] = pruneCertain(F_new,F_robust,F_xw,w);
for i = 1:length(F_new)
    if ~isempty(intersect(depends(F_new(i)),depends(w)))
        F_xw = F_xw + F_new(i);
    else
        F_robust = F_robust + F_new(i);
    end
end

function p = indexIn(x,y)
if ~isempty(x)
    for i = 1:length(x)
        p(i) = find(x(i)==y);
    end
else
    p = [];
end

function [F_x,F_w,F_xw,VariableType] = partitionModel(F,F_original,VariableType);

F_x = [];
F_w = [];
F_xw = [];

% x-var w_var aux_xw aux_w
if ~(isempty(VariableType.aux_with_w_dependence) & isempty(VariableType.aux_with_only_w_dependence))
    Dependency = spalloc(length(F_original),4,length(F));
    
    for i = 1:length(F_original)
        varF = depends(F_original(i));
        Dependency(i,1) = any(ismember(varF,VariableType.x_variables));
        Dependency(i,2) = any(ismember(varF,VariableType.w_variables));
        Dependency(i,3) = any(ismember(varF,VariableType.aux_with_w_dependence));
        Dependency(i,4) = any(ismember(varF,VariableType.aux_with_only_w_dependence));
    end
    
    LiftedUncertaintiesDescription = find(Dependency(:,1) == 0 & Dependency(:,3)==0);
    if ~isempty(LiftedUncertaintiesDescription)
        %     for i = LiftedUncertaintiesDescription(:)'
        %        vars = depends(F_original(i));
        %        vars = intersect(vars,VariableType.aux_with_only_w_dependence);
        %
        %     end
        reclassifyAsUncertain = depends(F_original(LiftedUncertaintiesDescription));
        [notused,reclassifyAsUncertain] = find(VariableType.Graph(reclassifyAsUncertain,:));
        reclassifyAsUncertain = unique(reclassifyAsUncertain);
        reclassifyAsUncertain = intersect(unique(reclassifyAsUncertain),getvariables(F));
        VariableType.aux_with_only_w_dependence = setdiff(VariableType.aux_with_only_w_dependence,reclassifyAsUncertain);
        VariableType.w_variables = union(VariableType.w_variables,reclassifyAsUncertain);
        
        VariableType.aux_with_w_dependence = union(VariableType.aux_with_w_dependence,VariableType.aux_with_only_w_dependence);
        VariableType.aux_with_only_w_dependence = [];
    end
end


% x-var w_var aux_xw aux_w
Dependency = spalloc(length(F),4,length(F));

for i = 1:length(F)
    varF = depends(F(i));
    Dependency(i,1) = any(ismember(varF,VariableType.x_variables));
    Dependency(i,2) = any(ismember(varF,VariableType.w_variables));
    Dependency(i,3) = any(ismember(varF,VariableType.aux_with_w_dependence));
    Dependency(i,4) = any(ismember(varF,VariableType.aux_with_only_w_dependence));
end

pureX = find(Dependency*[1;2;4;8] == 1);
pureW = find(Dependency(:,1) == 0 & Dependency(:,3)==0);
mixedXW = find(Dependency(:,1) | Dependency(:,3));
%mixedXW = find((Dependency(:,1) &  | Dependency(:,3));
mixedXW = setdiff(1:size(Dependency,1),union(pureW,pureX));

F_x = F(pureX);
F_w = F(pureW);
F_xw = F(mixedXW);

reclassifyAsUncertain = depends(F_w);
VariableType.aux_with_only_w_dependence = setdiff(VariableType.aux_with_only_w_dependence,reclassifyAsUncertain);
VariableType.w_variables = union(VariableType.w_variables,reclassifyAsUncertain);

function [VariableType,h_fixed,F_xw] = reformatObjective(h,F_xw,VariableType)
% Some pre-calc
x = recover(VariableType.x_variables);
w = recover(VariableType.w_variables);
xw = [x;w];
xind = find(ismembcYALMIP(getvariables(xw),getvariables(x)));
wind = find(ismembcYALMIP(getvariables(xw),getvariables(w)));
% Analyze the objective and try to rewrite any uncertainty into the format
% assumed by YALMIP
if ~isempty(h)
    
    [Q,c,f,dummy,nonquadratic] = vecquaddecomp(h,xw);
    Q = Q{1};
    c = c{1};
    f = f{1};
    
    if nonquadratic
        error('Objective can be at most quadratic, with the linear term uncertain');
    end
    
    Q_ww = Q(wind,wind);
    Q_xw = Q(xind,wind);
    Q_xx = Q(xind,xind);
    c_x = c(xind);
    c_w = c(wind);
    
    if nnz(Q_ww) > 0
        error('Objective can be at most quadratic, with the linear term uncertain');
    end
    % Separate certain and uncertain terms, place uncertain terms in the
    % constraints instead
    if is(h,'linear')
        if isempty(intersect(getvariables(w),getvariables(h)))
            h_fixed = h;
        else
            sdpvar t
            F_xw = F_xw + (h <= t);
            h_fixed = t;
            x = [x;t];
        end
    else
        h_fixed     = x'*Q_xx*x + c_x'*x + f;
        h_uncertain = 2*w'*Q_xw'*x + c_w'*w;
        if ~isa(h_uncertain,'double')
            sdpvar t
            F_xw = F_xw + (h_uncertain <= t);
            h_fixed = h_fixed + t;
            x = [x;t];
        end
    end
else
    h_fixed = [];
end
VariableType.x_variables = getvariables(x);