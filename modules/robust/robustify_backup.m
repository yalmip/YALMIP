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
% the robustification to be tractable. Please refer to the YALMIP Wiki for
% the current assumptions.
%
% Some options for the robustification strategies can be altered via the
% solver tag 'robust' in sdpsettings
%  'robust.lplp'  : Controls how linear constraints with affine
%                   parameterization in an uncertainty with polytopic
%                   description is handled. Can be either 'duality' or
%                   'enumeration'.
%  'polya'        : Controls the relaxation order of polynomials. If set to
%                   NAN, the polynomials will be eliminated by forcing the
%                   coefficients to zero
%
% See also UNCERTAIN

% Author Johan Löfberg
% $Id: robustify.m,v 1.55 2010-03-10 15:19:05 joloef Exp $

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

% Uncertainty in objective? If so, move to constraints by epigraph
if ~isempty(intersect(depends(h),depends(w)))
    t = sdpvar(1);
    F = [F,h < t];
    h = t;
end

% We keep track of auxilliary generated variables
nInitial = yalmip('nvars');

% Do we have any uncertain variables?
if isempty(w)
    unc_declarations = is(F,'uncertain');
    if any(unc_declarations)
        w = recover(getvariables(sdpvar(F(find(unc_declarations)))));
        F = F(find(~unc_declarations));
    else
        error('There is no uncertainty definition in the model.')
    end
end

% Figure out which variables are uncertain, certain, and lifted variables
% in the uncertainty description (this code is slow and buggy as ....)
%[x,w,x_variables,w_variables,aux_variables,F,failure] = robust_classify_variables_new(F,h,ops,w);

[x,w,x_variables,w_variables,aux_variables,F,failure,VariableClassifications] = robust_classify_variables(F,h,ops,w);
if failure
    return
end

if ops.verbose
    disp('***** Starting YALMIP robustification module. *****');
    disp([' - Detected ' num2str(length(w)) ' uncertain variables']);
end

% Integer variables are OK in x, but not in the uncertainty (robustification
% is based on strong duality in w-space)
integervars = [yalmip('binvariables') yalmip('intvariables')];
ind = find(is(F,'integer') | is(F,'binary'));
if ~isempty(ind)
    integervars = [integervars getvariables(F(ind))];
    if any(ismember(w_variables,integervars))
        failure = 1;
        return
    end
end

% Find  uncertainty description, uncertain and certain constraints
F_w = set([]);
F_x = set([]);
F_xw = set([]);
% FIXME: Is this really correct
x2_variables = setdiff(x_variables,yalmip('auxvariablesW'));

for i = 1:length(F)
    if all(ismember(depends(F(i)),w_variables))
        % Uncertainty definition
        F_w = F_w + F(i);
    elseif all(ismember(depends(F(i)),x2_variables))
        % Certain constraint
        F_x = F_x +  F(i);
    else
        % Uncertain constraint
        F_xw = F_xw + F(i);
    end
end

% Limitation in the modelling language...
if ~isempty(intersect(intersect(depends(F_xw),depends(F_w)),aux_variables))
   disp('You are most likely using a mixture of operator to describe the');
   disp('uncertainty set. This is currently not supported.');
   disp('Please model the constraint manually.');
   error('Uncertain model does not satisfy assumptions (nonlinear operator on uncertainty in uncertain constraint)');
end

if length(F_w)==0
    error('There is no uncertainty description in the model.');
end

% Some pre-calc
xw = [x;w];
xind = find(ismembc(getvariables(xw),getvariables(x)));
wind = find(ismembc(getvariables(xw),getvariables(w)));
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
            F_xw = F_xw + set(h < t);
            h_fixed = t;
            x = [x;t];
        end
    else
        h_fixed     = x'*Q_xx*x + c_x'*x + f;
        h_uncertain = 2*w'*Q_xw'*x + c_w'*w;
        if ~isa(h_uncertain,'double')
            sdpvar t
            F_xw = F_xw + set(h_uncertain < t);
            h_fixed = h_fixed + t;
            x = [x;t];
        end
    end
else
    h_fixed = [];
end

% Convert quadratic constraints in uncertainty model to SOCPs. This will
% enable us to use duality based removal of uncertainties in linear
% inequalities. We keep information about quadratic expression though, in
% order to use them if possible in 2-norm explicit maximization
F_w = convertquadratics(F_w);

% Export uncertainty model to numerical format
ops.solver = '';
[aux1,aux2,aux3,Zmodel] = export(F_w,[],ops,[],[],1);

if ~isempty(Zmodel)
    if length(Zmodel.c) ~= length(w)
        error('Some uncertain variables are unconstrained.')
    end
else
    error('Failed when exporting a model of the uncertainty.')
end

% The uncertainty model is in the full w-space. However, it might
% happen that the uncertainty is separable. Find groups of uncertain
% variables that can/should be treated separately.
uncertaintyGroups = findSeparate(Zmodel);

if ops.verbose
    disp([' - Detected ' num2str(length(uncertaintyGroups)) ' independent group(s) of uncertain variables']);
end

% ************* EXPERIMENTAL **********************************************
% Trying to take care of cases such as norm([x+w;x-w]), i.e. epigraphs with
% uncertainty.
if intersect(yalmip('auxvariablesW'),x_variables)
    disp([' - Detected ' num2str(length(intersect(x_variables,yalmip('auxvariablesW')))) ' uncertainty dependent auxilliary variables']);
    if isequal(Zmodel.K.s,0) & isequal(Zmodel.K.q,0)  &  ~isempty(strfind('enumeration',ops.robust.auxreduce)) %isequal(ops.robust.auxreduce,'enumeration')
        disp([' - Using enumeration due to the uncertainty dependent auxilliary variables']);
        ops.robust.lplp = 'enumeration';
        forced_enumeration = 1;
    elseif ~(any(is(F_xw,'lmi')) | any(is(F_xw,'socp'))) | ~isempty(strfind('projection',ops.robust.auxreduce))
        if isequal(ops.robust.auxreduce,'projection')
            disp([' - Projecting out auxilliary variables. This can take time...']);
        else
            disp([' - Forced to project out auxilliary variables due to conic uncertainty. This can take time...']);
        end
        F_xw = projectOut(F_xw,w,unique(yalmip('auxvariablesW')),uncertaintyGroups,ops);
        forced_enumeration = 0;
    else
        disp([' - Cannot robustify model with uncertainty dependent auxilliary variables and conic uncertainty']);
        error('robustification failed');
    end
else
    forced_enumeration = 0;
end

% *************************************************************************


% Separate the uncertainty models accroding to uncertainty groups
separatedZmodel = separateUncertainty(Zmodel,uncertaintyGroups);

% Conversion of bounded variables that have been modelled using
% the norm operator (remove the epigraph variable to ensure explicit
% maximization is used). This will be generalized in the next version
[separatedZmodel, uncertaintyGroups] = convertUncertaintyGroups(separatedZmodel,uncertaintyGroups);

% Code will be added to detect uncertainty cases in a more general and
% modular way. Additional code will also be added to find hidden simple
% structures, such as norm(w,1)<1, which currently is treated as a general
% polytopic uncertainty, since the expansion hides the simplicity
% 'Bounds', 'Simplex', 'Conic', 'Polytopic', '2-norm', '1-norm', 'inf-norm'
[uncertaintyTypes,separatedZmodel,uncertaintyGroups] = classifyUncertainty(separatedZmodel,uncertaintyGroups);

% OK, we are done with the initial analysis of the involved variables, and
% check of the objective function.

% At this point, we have to decide on the algorithm we should use for
% robustifying the constraints. There are a couple of alternatives,
% depending on uncertainty and constraints
% 1. Polya:       Polynomial uncertainty dependence, simplex uncertainty,
%                 can only be applied on LP constraints
% 2. Elimination: Last resort, tries to cancel all nonlinear uncertainties
%                 by setting coefficients to zero
% 3. Explicit:    Linear uncertainty dependence, box-model uncertainty, can
%                 only be applied on LP constraints
% 4a. Enumeration:Linear uncertainty dependence, polytopic uncertainty,
%                 arbitrary type of constraints (convex)
% 4b. Duality:    Linear uncertainty dependence, conic uncertainty, can
%                 only be applied on LP constraints

% Initialize a set of constraints that correspond to robustified
% constraints, and constraints related to Polyas theorem, and elimination
% of nonlinear terms and uncertain equalities
F_robust = set([]);

% We begin by checking to see if the user wants to apply Polyas theorem.
% If that is the case, search for simplex structures, and apply Polyas.
if ~isnan(ops.robust.polya) & any(strcmp(uncertaintyTypes,'simplex')) & ~forced_enumeration
    F_polya = [];
    % Recursively apply Polya relaxation w.r.t each simplex
    for i = find(strcmp(uncertaintyTypes,'simplex'))
        [F_xw, F_polya] = filter_polya(F_xw+F_polya,w(uncertaintyGroups{i}),ops.robust.polya);
    end
    [F_xw,F_robust] = pruneCertain(F_polya,F_robust,F_xw,w);
end

% LP constraints with quadratic dependence and quadratic uncertainty region
% can be handled tightly using the S-procedure
if all(strcmp(uncertaintyTypes,'2-norm')) & length(uncertaintyTypes)==1 & ~forced_enumeration
    [F_xw,F_sprocedure] = filter_sprocedure(F_xw,w,separatedZmodel,ops);
    F_robust = F_robust + F_sprocedure;
end

% There might still be nonlinearities left in the model. These have to be
% removed. We remove all terms with w-degree larger than 1
[F_xw,F_elimination] = filter_eliminatation(F_xw,w,1);
F_robust = F_robust + F_elimination;

% Equality constraints cannot be part of an uncertain problem. Any
% dependence w.r.t w in inequalities has to be removed
F_eq = F_xw(find(is(F_xw,'equality')));
F_xw = F_xw - F_eq;
[F_eq_left,F_eliminate_equality] = filter_eliminatation(F_eq,w,0);
F_robust = F_robust + F_eliminate_equality + F_eq_left;

% The problem should now be linear in the uncertainty, with no uncertain
% equality constraints. Hence, now we apply explicit maximization,
% enumeration or duality-based robustification.

% We start with the norm balls
if ~forced_enumeration
    for i = 1:length(uncertaintyTypes)
        %if strcmp(uncertaintyTypes{i},'1-norm') | strcmp(uncertaintyTypes{i},'2-norm') | strcmp(uncertaintyTypes{i},'inf-norm')
        if ismember(uncertaintyTypes{i},{'1-norm','2-norm','inf-norm'})
            F_lp = F_xw(find(is(F_xw,'elementwise')));
            F_xw = F_xw - F_lp;
            F_flt = filter_normball(F_lp,separatedZmodel{i},x,w(uncertaintyGroups{i}),w,uncertaintyTypes{i},ops);
            [F_xw,F_robust] = pruneCertain(F_flt,F_robust,F_xw,w);
        end
    end
end

% Pick out the uncertain linear equalities and robustify using duality if
% user has opted for this or the uncertainty is conic.
conic = ~isequal(Zmodel.K.s,0) | ~isequal(Zmodel.K.q,0);
if (conic | isequal(ops.robust.lplp,'duality')) & ~forced_enumeration
    F_lp = F_xw(find(is(F_xw,'elementwise')));
    F_xw = F_xw - F_lp;
    F_robust = F_robust + filter_duality(F_lp,Zmodel,x,w);
end

% Robustify remaining uncertain LP/SOCP/SDP constraints and robustify by
% enumeration.
F_conic = F_xw(find(is(F_xw,'sdp') | is(F_xw,'socc') | is(F_xw,'elementwise')));
F_xw = F_xw - F_conic;
[F_temp,enumerationfailed] = filter_enumeration(F_conic,Zmodel,x,w,ops,uncertaintyTypes,separatedZmodel);
if enumerationfailed
    % Reset to previous state
    F_xw = F_xw + F_conic;
else
    F_robust = F_robust + F_temp;
end

if enumerationfailed
    % Enumeration failed, probably due to lack of MPT. If problem is conic,
    % we are in troubles. If simple LP, we can use duality
    if conic
        disp(' - Enumeration of uncertainty polytope failed. Missing Multiparametric Toolbox?')
        error('Enumeration failed (lacking MPT?),  and due to conic constraints, duality cannot be used');
    else
        F_lp = F_xw(find(is(F_xw,'elementwise')));
        F_xw = F_xw - F_lp;
        nv = yalmip('nvars');
        F_filter = filter_duality(F_lp,Zmodel,x,w);
        disp([' - Duality introduced ' num2str(yalmip('nvars')-nv') ' variables, ' num2str(nnz(is(F_filter,'equality'))) ' equalities, ' num2str(nnz(is(F_filter,'elementwise'))) ' LP inqualities and ' num2str(nnz(is(F_filter,'sdp'))+nnz(is(F_filter,'socp'))) ' conic constraints']);
        F_robust = F_robust + F_filter;
    end
end

% If there is anything left now, it means that we do not support it (such
% as conic uncertainty in conic constraint)
if length(F_xw) > 0
    if any(~islinear(F_xw))
        error('There are some uncertain constraints which cannot be robustified by YALMIP')
    else
        F_robust = F_robust + F_xw;
    end
end

% Return the robustfied model
F = F_robust+F_x;
h = h_fixed;

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



function groups = findSeparate(model)

% This is an early bail-out to avoid any errors during the devlopment of
% new features. Separable constraints are only support for Polya models
if any(model.K.s > 0) %| any(model.K.q >0)
    groups{1} = 1:size(model.F_struc,2)-1;
    return
end

X = zeros(size(model.F_struc,2)-1);
top  = 1;
if model.K.f + model.K.l > 0
    A = model.F_struc(top:model.K.f+model.K.l,2:end);
    for i = 1:size(A,1)
        X(find(A(i,:)),find(A(i,:))) = 1;
    end
    top = top + model.K.f + model.K.l;
end

if any(model.K.q > 0)
    for j = 1:length(model.K.q)
        A = model.F_struc(top:top+model.K.q(j)-1,2:end);top = top + model.K.q(j);
        A = sum(abs(A),1);
        for i = 1:size(A,1)
            X(find(A),find(A)) = 1;
        end
    end
end

if any(model.K.s > 0)
    for j = 1:length(model.K.s)
        A = model.F_struc(top:top+model.K.s(j)^2-1,2:end);top = top + model.K.s(j)^2;
        A = sum(abs(A),1);
        for i = 1:size(A,1)
            X(find(A),find(A)) = 1;
        end
    end
end

[a,b,c,d] = dmperm(X);
for i = 1:length(d)-1
    groups{i} = sort(a(d(i):d(i+1)-1));
end

function Zmodel = convertuncertainty(Zmodel);
% Temoporary hack, will be generalized once the framework for multiple
% uncertainty models is supported
% We are looking for k>t, -tw<t
if size(Zmodel,1) == 1+(size(Zmodel,2)-1)*2 & Zmodel.K.f==0 & Zmodel.K.l==size(Zmodel.F_struc,1)
    n = size(Zmodel.F_struc,2)-1;
    if isequal(Zmodel.F_struc(:,2:end),sparse([zeros(1,n-1) -1;[eye(n-1);-eye(n-1)] ones(2*(n-1),1)]))
        Zmodel.F_struc = [ones(2*n,1)*Zmodel.F_struc(1,1) [eye(n);-eye(n)]];
        Zmodel.K.l = 2*n;
    end
end


function [separatedZmodel, uncertaintyGroups] = convertUncertaintyGroups(separatedZmodel,uncertaintyGroups);
% Temoporary hack, will be generalized once the framework for multiple
% uncertainty models is supported

% We are looking for k>t, -t<w<t. This is the slightly redundant model for
for i = 1:length(separatedZmodel)
    % -k<w<k genrated when YALMIP encounters norm(w,inf)<k
    if size(separatedZmodel{i},1) == 1+(size(separatedZmodel{i},2)-1)*2 & separatedZmodel{i}.K.f==0 & separatedZmodel{i}.K.l==size(separatedZmodel{i}.F_struc,1)
        n = size(separatedZmodel{i}.F_struc,2)-1;
        if isequal(separatedZmodel{i}.F_struc(:,2:end),sparse([zeros(1,n-1) -1;[eye(n-1);-eye(n-1)] ones(2*(n-1),1)]))
            k = separatedZmodel{i}.F_struc(1,1);
            c = separatedZmodel{i}.F_struc(2:2+n-2,1);
            separatedZmodel{i}.F_struc = [[k-c;k;c+k;k] [-eye(n);eye(n)]];
            separatedZmodel{i}.K.l = 2*n;
        end
        
    elseif separatedZmodel{i}.K.l == 1 & separatedZmodel{i}.K.q(1)>0 & length(separatedZmodel{i}.K.q(1))==1
        % Could be norm(w,2) < r.
        [n,m] = size(separatedZmodel{i}.F_struc);
        if n==m & n>=3
            n = n-2;
            
            if isequal(separatedZmodel{i}.F_struc(:,2:end),sparse([zeros(2,n) [-1;1];eye(n) zeros(n,1)]))% isequal(separatedZmodel{i}.F_struc(:,2:end),sparse([zeros(2,n) [-1;1];eye(n) zeros(n,1)]))
                if separatedZmodel{i}.F_struc(1,1)>0 %& nnz(separatedZmodel{i}.F_struc(2:end,1))==0
                    % the user has written norm(w-xc) < r. YALMIP will handle
                    % this using the nonlinear operator framework, and
                    % introduce a variable t and write it as norm(w)<t, t<r
                    % The variable t will be appended to the list of
                    % uncertain variables, and thus the model is a SOCP
                    % uncertainty in (w,t). We however want to interpret it
                    % as a simple quadratic model. We thus write it as
                    % (w,t)Q(w,t) < 1. Since t actually is unbounded in
                    % this form (since t is an auxilliary variable that we
                    % now through away, we arbitrarily say Q=[I/r^2 0;0 1]
                    r = separatedZmodel{i}.F_struc(1,1);
                    center = -separatedZmodel{i}.F_struc(3:end,1);
                    % separatedZmodel{i}.Q = blkdiag(r*eye(n),1);
                    % separatedZmodel{i}.center = [center;0];
                    
                    separatedZmodel{i}.K.l = 0;
                    separatedZmodel{i}.K.q = n+1;
                    % r = separatedZmodel{i}.F_struc(1,1);
                    % center = -separatedZmodel{i}.F_struc(3:end,1);
                    separatedZmodel{i}.F_struc = [[1+r^2;-2*center;1-r^2] [zeros(1,n);eye(n)*2;zeros(1,n)]];
                    uncertaintyGroups{1}(end)=[];
                end
            end
        end
    end
    
end



% Separate the uncertainty models
function newModels = separateUncertainty(Zmodel,uncertaintyGroups);

for i = 1:length(uncertaintyGroups)
    data = Zmodel.F_struc(:,[1 1+uncertaintyGroups{i}]);
    
    data_f = data(1:Zmodel.K.f,:);
    data_l = data(Zmodel.K.f+1:Zmodel.K.f+Zmodel.K.l,:);
    %data_q = data(Zmodel.K.f+Zmodel.K.l+1:sum(Zmodel.K.q),:);
    %data_s = data(Zmodel.K.f+Zmodel.K.l+1:sum(Zmodel.K.q)+1:end,:);
    
    eqIndex = find(any(data_f(:,2:end),2));
    liIndex = find(any(data_l(:,2:end),2));
    
    newModels{i}.F_struc = [data_f(eqIndex,:);data_l(liIndex,:)];
    newModels{i}.K.f = length(eqIndex);
    newModels{i}.K.l = length(liIndex);
    newModels{i}.K.q = [];
    top = Zmodel.K.f+Zmodel.K.l+1;
    if Zmodel.K.q(1)>0
        for j = 1:length(Zmodel.K.q)
            data_q = data(top:top+Zmodel.K.q(j)-1,:);
            if nnz(data_q(:,2:end))>0
                newModels{i}.K.q(end+1) = Zmodel.K.q(j);
                newModels{i}.F_struc = [newModels{i}.F_struc;data_q];
            end
            top = top + Zmodel.K.q(j);
        end
    end
    if isempty(newModels{i}.K.q)
        newModels{i}.K.q = 0;
    end
    newModels{i}.K.s = [];
    top = Zmodel.K.q+Zmodel.K.f+Zmodel.K.l+1;
    if Zmodel.K.s(1)>0
        for j = 1:length(Zmodel.K.s)
            data_s = data(top:top+Zmodel.K.s(j)^2-1,:);
            if nnz(data_s(:,2:end))>0
                newModels{i}.K.s(end+1) = Zmodel.K.s(j);
                newModels{i}.F_struc = [newModels{i}.F_struc;data_s];
            end
            top = top + Zmodel.K.q(j);
        end
    end
    if isempty(newModels{i}.K.s)
        newModels{i}.K.s = 0;
    end
    
    newModels{i}.K.s = 0;
    newModels{i}.variabletype = Zmodel.variabletype(uncertaintyGroups{i});
end


function [uncertaintyTypes,separatedZmodel,uncertaintyGroups] = classifyUncertainty(separatedZmodel,uncertaintyGroups)

for i = 1:length(separatedZmodel)
    uncertaintyTypes{i} = 'unclassified';
end

% Look for simplicies, which can be used in Polya
simplex_members = find_simplex_models(separatedZmodel);
for i = find(simplex_members)
    uncertaintyTypes{i} = 'simplex';
end

% Look for simple bounds, and them combine them
for i = 1:length(separatedZmodel)
    if strcmp(uncertaintyTypes{i},'unclassified')
        if any(separatedZmodel{i}.K.q > 0) | any(separatedZmodel{i}.K.s > 0)
            simplebounds = 0;
        else
            [aux,lower,upper] = find_simple_variable_bounds(separatedZmodel{i});
            simplebounds = ~isinf(lower) & ~isinf(upper);
        end
        if all(simplebounds)
            if aux.K.l == 0
                uncertaintyTypes{i} = 'inf-norm';
                separatedZmodel{i}.lb = lower;
                separatedZmodel{i}.ub = upper;
            end
        end
    end
end

j = find(strcmp(uncertaintyTypes,'inf-norm'));
if length(j)>1
    allBounded = [];
    lb = [];
    ub = [];
    for i = 1:length(j)
        allBounded = [allBounded; uncertaintyGroups{j(i)}];
        lb = [lb;separatedZmodel{j(i)}.lb];
        ub = [ub;separatedZmodel{j(i)}.ub];
        if any(lb > ub)
            error('There are inconsistent bounds in the uncertainty model');
        end
    end
    separatedZmodel{j(1)}.lb = lb;
    separatedZmodel{j(1)}.ub = ub;
    uncertaintyGroups{j(1)} = allBounded(:)';
    keep = setdiff(1:length(separatedZmodel),j(2:end));
    separatedZmodel = {separatedZmodel{keep}};
    uncertaintyGroups = {uncertaintyGroups{keep}};
    uncertaintyTypes = {uncertaintyTypes{keep}};
end

% Look for 2-norm balls norm(w) < r, written using norm(w) or w'*w
for i = 1:length(separatedZmodel)
    if strcmp(uncertaintyTypes{i},'unclassified')
        if ~any(separatedZmodel{i}.K.s > 0) &  separatedZmodel{i}.K.l==0 & separatedZmodel{i}.K.f==0 & length(separatedZmodel{i}.K.q)==1
            % Hmm, only 1 SOCP
            if all(separatedZmodel{i}.F_struc(1,2:end)==0) % r > norm(f)
                if isequal(separatedZmodel{i}.F_struc(2:end,2:end),speye(separatedZmodel{i}.K.q(1)-1))
                    uncertaintyTypes{i} = '2-norm';
                    separatedZmodel{i}.center = -separatedZmodel{i}.F_struc(2:end,1);
                    separatedZmodel{i}.r = sqrt(max(0,separatedZmodel{i}.F_struc(1,1)^2));
                elseif size(separatedZmodel{i}.F_struc,1)==size(separatedZmodel{i}.F_struc,2)
                    S = separatedZmodel{i}.F_struc(2:end,2:end);
                    [ii,jj,kk] = find(S);
                    if all(kk==2) & length(ii)==length(unique(ii)) & length(jj)==length(unique(jj))
                        uncertaintyTypes{i} = '2-norm';
                        separatedZmodel{i}.center = -separatedZmodel{i}.F_struc(2:end,1)/2;
                        separatedZmodel{i}.r = sqrt((separatedZmodel{i}.F_struc(1,1)^2-separatedZmodel{i}.F_struc(end,1)^2)/4);
                    end
                elseif separatedZmodel{i}.F_struc(1,1)^2-separatedZmodel{i}.F_struc(end,1)^2>0
                    S = separatedZmodel{i}.F_struc(2:end,2:end);
                    [ii,jj,kk] = find(S);
                    if all(kk==2) & length(ii)==length(unique(ii)) & length(jj)==length(unique(jj))
                        uncertaintyTypes{i} = '2-norm';
                        separatedZmodel{i}.center = -separatedZmodel{i}.F_struc(2:end-1,1)/2;
                        separatedZmodel{i}.r = sqrt((separatedZmodel{i}.F_struc(1,1)^2-separatedZmodel{i}.F_struc(end,1)^2)/4);
                    end
                end
            end
        end
    end
end

% Look for 1-norm balls |w|1 < r^2
% We are looking for the auto-generate model -t<w<t, sum(t)<s,s<r^2
for i = 1:length(separatedZmodel)
    if strcmp(uncertaintyTypes{i},'unclassified')
        if ~any(separatedZmodel{i}.K.s > 0) &  separatedZmodel{i}.K.l>0 & separatedZmodel{i}.K.f==0 & ~any(separatedZmodel{i}.K.q > 0)
            %if all(separatedZmodel{i}.F_struc(2:end,1)==0)
            if separatedZmodel{i}.F_struc(1,1)>0
                n = (size(separatedZmodel{i}.F_struc,2)-2)/2;
                if n==fix(n)
                    try
                        if all(separatedZmodel{i}.F_struc(n+2:end-1,1)==-separatedZmodel{i}.F_struc(2:2+n-1,1))
                            shouldbe = [zeros(1,n) -1 zeros(1,n);eye(n) zeros(n,1) eye(n);-eye(n) zeros(n,1) eye(n);zeros(1,n) 1 -ones(1,n)];%; zeros(n,n+1)  eye(n)];;zeros(1,n) 1 zeros(1,n)];
                            if  isequal(full(separatedZmodel{i}.F_struc(:,2:end)),shouldbe)
                                %if all(separatedZmodel{i}.F_struc(2:end,1)==0)
                                uncertaintyTypes{i} = '1-norm';
                                separatedZmodel{i}.r = separatedZmodel{i}.F_struc(1,1);
                                separatedZmodel{i}.center =separatedZmodel{i}.F_struc(n+2:end-1,1);
                                %end
                            end
                        end
                    catch
                        %FIX ME, caught by -1<w<1, sum(w)<1
                    end
                end
            end
            %end
        end
    end
end


% Look for polytopes
for i = 1:length(separatedZmodel)
    if strcmp(uncertaintyTypes{i},'unclassified')
        if ~any(separatedZmodel{i}.K.q > 0) | any(separatedZmodel{i}.K.s > 0)
        end
    end
end



% Look for conic representations
function [F_xw,F_robust] = pruneCertain(F_new,F_robust,F_xw,w);
for i = 1:length(F_new)
    if ~isempty(intersect(depends(F_new(i)),depends(w)))
        F_xw = F_xw + F_new(i);
    else
        F_robust = F_robust + F_new(i);
    end
end





function Fnew = projectOut(F,w,newAuxVariables,uncertaintyGroups,ops)

w_variables = getvariables(w);

if 1
    Graph = yalmip('getdependence');
    for i = 1:length(uncertaintyGroups)
        Graph(w_variables(uncertaintyGroups{i}),w_variables(uncertaintyGroups{i}))=1;
    end
    aux_and_w = union(newAuxVariables,w_variables);
    F_lp = F(find(is(F,'elementwise')));
    X = sdpvar(F_lp);
    Xvars = getvariables(X);
    G = getbase(X);G = G(:,2:end);
    for i = 1:size(G,1)
        j = Xvars(find(G(i,:)));
        Graph(j,j) = 1;
    end
    Graph(:,setdiff(1:size(Graph,2),aux_and_w))=0;
    Graph(setdiff(1:size(Graph,2),aux_and_w),:)=0;
    
    
 if 0   
    [ss,cc] = graphconncomp(sparse(Graph),'DIRECTED',false);
    cc = cc(newAuxVariables);
    uniqueGroups = unique(cc);
    for r = 1:length(uniqueGroups)
        commonProjections{r} = newAuxVariables(find(cc == uniqueGroups(r)));
    end
 else
     Graph2 = Graph(1:max(aux_and_w),1:max(aux_and_w));
     Graph2 = Graph2 + speye(length(Graph2));
     [aa,bb,cc,dd] = dmperm(Graph2);
     commonProjections = {};
     for r = 1:length(dd)-1
         comps = dd(r):dd(r+1)-1;
         comps = intersect(aa(comps),newAuxVariables);
         if ~isempty(comps)
            commonProjections{end+1}  = comps;
         end         
     end
 end
    
else
    % Build a graph showing dependence of auxilliary variables depending on
    % uncertainty, and uncertain variables. for instance -s < x+w < s, hence s
    % depends on w
    
    
    Graph = sparse(0);
    for i = 1:length(F)
        F_vars = getvariables(F(i));
        var_w = intersect(F_vars,w_variables);
        var_aux = intersect(F_vars,newAuxVariables);
        Graph(indexIn(var_w,w_variables),indexIn(var_aux,newAuxVariables)) = 1;
    end
    
    
    
    % alls = union(w_variables,newAuxVariables);
    % Graph = Graph(newAuxVariables,:);
    % for i = 1:size(Graph,1)
    %     j = ismember(Graph(i,:),newAuxVariables)
    %     while ~isempty(j)
    %         for k = j(:)'
    %
    %         end
    %         any(Graph(i,:)
    % Graph = Graph(newAuxVariables,union(w_variables,newAuxVariables))
    % for i = newAuxVariables
    %     Graph(i,:)
    % end
    
    % we might have things like s1+s2 < t1 where t1 is an epigrapg variable.
    % That means t1 depends on w too
    for i = 1:length(F)
        F_vars = getvariables(F(i));
        var_w = intersect(F_vars,w_variables);
        if isempty(var_w) & all(ismember(F_vars,newAuxVariables))
            var_aux = intersect(F_vars,newAuxVariables);
            k = indexIn(var_aux,newAuxVariables);
        end
        Graph(indexIn(var_w,w_variables),indexIn(var_aux,newAuxVariables)) = 1;
    end
    
    ExtGraph = zeros(length(w_variables) + length(newAuxVariables));
    for i = 1:length(uncertaintyGroups)
        ExtGraph(uncertaintyGroups{i},uncertaintyGroups{i})=1;
    end
    
    for j = 1:size(Graph,1)
        ExtGraph(j,length(w_variables)+find(Graph(j,:))) = 1;
        ExtGraph(length(w_variables)+find(Graph(j,:)),j) = 1;
    end
    
    [ss,cc] = graphconncomp(sparse(ExtGraph),'DIRECTED',false);
    cc = cc(length(w_variables)+1:end);
    cc_red = cc(find(any(full(Graph),1)));
    
    for r = 1:length(unique(cc_red))
        commonProjections{r} = newAuxVariables(find(cc == cc_red(r)));
    end
    
    commonProjections = {};
    for i = 1:length(uncertaintyGroups)
        [ii,jj] = find(Graph(uncertaintyGroups{i},:));
        if ~isempty(jj)
            commonProjections{end+1} = newAuxVariables(unique(jj));
        end
    end
    
end
if ops.verbose
    disp([' - * Detected ' num2str(nnz(any(Graph,1))) ' graph variables depending on uncertainty']);
    disp([' - * Partitioned these into ' num2str(length(commonProjections)) ' group(s)']);
end

keep = ones(length(F),1);
Fnew = [];
started = 0;
for clique = 1:length(commonProjections)
    F_lp = [];
    for i = 1:length(F)
        F_vars = getvariables(F(i));
        var_w = intersect(F_vars,w_variables);
        var_aux = intersect(F_vars,commonProjections{clique});
        
        if is(F(i),'elementwise')
            if any(ismember(getvariables(F(i)),commonProjections{clique}))
                F_lp = [F_lp, F(i)];
                keep(i) = 0;
                if ~started
                    started = 1;
                    if ops.verbose
                        disp([' - * Performing projections of uncertain graph variables...']);
                    end
                end
            end
        end
    end
    if ~isempty(F_lp)
        X = sdpvar(F_lp);
        [Ab] = getbase(X);
        b = Ab(:,1);
        A = -Ab(:,2:end);
        allVariables = getvariables(X);
        newAuxVariables = intersect(commonProjections{clique},allVariables);
        j = [];
        for i = 1:length(newAuxVariables)
            j(i) = find(allVariables==newAuxVariables(i));
        end
        
        [A,b] = fourierMotz(A,b,sort(j));
        A(:,j)=[];
        left = recover(setdiff(allVariables,newAuxVariables));
        X = b-A*left(:);
        Fnew = [Fnew,[[X > 0] : 'Projected uncertain']];
    end
end
if started & ops.verbose
    disp([' - * Done with projections (generated ' num2str(length(sdpvar(Fnew))) ' constraints)']);
end
Fnew = [Fnew,F(find(keep))];

function [A,b] = fourierMotz(A,b,dim)
for i=dim
    positive=find(A(:,i) > 0);
    negative=find(A(:,i) < 0);
    null=find(A(:,i) == 0);
    
    nr=length(null) + length(positive)*length(negative);
    nc=size(A,1);
    C=sparse(zeros(nr,nc));
    
    AA=A(:,i);
    row=1;
    for j=(positive)'
        for k=(negative)'
            C(row,j)=-AA(k);
            C(row,k)=AA(j);
            row=row+1;
        end
    end
    
    for j=(null)'
        C(row,j)=1;
        row=row+1;
    end
    
    % Compute new Matrix P.H and P.K
    A=C*A;
    b=C*b;
end

function p = indexIn(x,y)
if ~isempty(x)
    for i = 1:length(x)
        p(i) = find(x(i)==y);
    end
else
    p = [];
end



