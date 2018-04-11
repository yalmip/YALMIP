function [sol,info] = solvebilevel(OuterConstraints,OuterObjective,InnerConstraints,InnerObjective,InnerVariables,options)
%SOLVEBILEVEL Simple global bilevel solver
%
%   min        OO(x,y)
%   subject to CO(x,y)>=0
%              y = arg min OI(x,y)
%              subject to CI(x,y)>=0
%
%   [DIAGNOSTIC,INFO] = SOLVEBILEVEL(CO, OO, CI, OI, y, options)
%
%   diagnostic : Struct with standard YALMIP diagnostics
%   info       : Bilevel solver specific information
%
%   Input
%      CO       : Outer constraints (linear elementwise)
%      OO       : Outer objective (convex quadratic)
%      CI       : Inner constraints (linear elementwise)
%      OI       : Inner objective (convex quadratic)
%      y        : Inner variables
%      options  : solver options from SDPSETTINGS.
%
%   The behaviour of the bilevel solver can be controlled
%   using the field 'bilevel' in SDPSETTINGS
%
%      bilevel.outersolver : Solver for outer problems with inner KKT removed
%      bilevel.innersolver : Solver for inner problem
%      bilevel.rootcut     : Number of cuts (based on complementary
%                            constraints) added in root (experimental)
%      bilevel.relgaptol   : Termination tolerance
%      bilevel.compslacktol: Tolerance for accepting complementary slackness
%      bilevel.feastol     : Tolerance for feasibility in outer problem
%
%
%   See also SDPVAR, SDPSETTINGS, SOLVESDP

% min f(x,y) s.t g(x,y)<0, y = argmin  [x;y]'*H*[x;y]+e'[x;y]+f, E[x;y]<d

if nargin<6
    options = sdpsettings;
elseif isempty(options)
    options = sdpsettings;
end

y = InnerVariables;

if ~isempty(InnerConstraints)
    if any(is(InnerConstraints,'sos2'))
        error('SOS2 structures not allowed in inner problem');
    end
end

% User wants to use fmincon, cplex or something like
if strcmp(options.bilevel.algorithm,'external')
    % Derive KKT conditions of inner problem, append with outer, and solve
    % using standard solver
    z = [depends(OuterConstraints) depends(OuterObjective) depends(InnerObjective) depends(InnerConstraints)];
    z = setdiff(z,depends(y));
    z = recover(unique(z));
    [K,details] = kkt(InnerConstraints,InnerObjective,z,options);    
    Constraints = [K,OuterConstraints];   
    options.solver = options.bilevel.outersolver;
    sol = solvesdp(Constraints,OuterObjective,options);
    info = [];
    return
end

% Export the inner model, and select solver
options.solver = options.bilevel.innersolver;
if isa(InnerObjective, 'double') || is(InnerObjective,'linear')
    [Imodel,Iax1,Iax2,inner_p] = export(InnerConstraints,InnerObjective,options,[],[],0);
elseif is(InnerObjective,'quadratic')
    % We have to be a bit careful about cases such as x'y. This is convex in
    % the inner problem, since x is constant there.
    % [Q,c,f,dummy,nonquadratic] = vecquaddecomp(InnerObjective);
    % Extract model for a fake quadratic model
    % [InnerConstraints,failure] = expandmodel(InnerConstraints,InnerObjective)
    %[Imodel,Iax1,Iax2,inner_p] = export(InnerConstraints,dummy'*diag(1+diag(Q{1}))*dummy+c{1}'*dummy,options,[],[],0);
   % toptions = options;
   % toptions.expandbilinear = 1;
    yy = recover(setdiff(depends(y),setdiff(depends(InnerObjective),depends(y))));
    [Imodel,Iax1,Iax2,inner_p] = export(InnerConstraints,yy'*yy+sum(recover(depends(InnerObjective))),options,[],[],0);
    [Q,c,f,dummy,nonquadratic] = vecquaddecomp(InnerObjective,recover(inner_p.used_variables));
    %[Imodel,Iax1,Iax2,inner_p] = export(InnerConstraints,InnerObjective,options,[],[],0);
    % Now plug in the real quadratic function
    if ~isequal(getvariables(dummy),inner_p.used_variables)
        error('This quadratic form is not supported yet. Please make feature request')
    else
        inner_p.Q = Q{1};
        inner_p.c = c{1};
        inner_p.f = f{1};
    end
else
    error('Only LPs or convex QPs allowed as inner problem');
end

% Modeling of inner problem might have lead to more decision variables in
% the inner problem. Append these
v1 = getvariables(y);
v2 = inner_p.used_variables(inner_p.extended_variables);
v3 = inner_p.used_variables(inner_p.aux_variables(:));
y = recover(unique([v1(:);v2(:);v3(:)]));

% Export the outer model, and select solver
options.solver = options.bilevel.outersolver;
options.bmibnb.diagonalize = 0;
[Omodel,Oax1,Oax2,outer_p] = export(OuterConstraints,OuterObjective,options,[],[],0);
if isstruct(Oax2)
   sol = Oax2;
   info = 2;
   return
end
% Export a joint model with KKT removed, to simplify some setup later
% [Auxmodel,Auxax1,Auxax2,outerinner_p] = export([OuterConstraints,InnerConstraints],OuterObjective+pi*InnerObjective,options,[],[],0);

if ~all(inner_p.variabletype==0) | ~isequal(inner_p.K.s,0) | ~isequal(inner_p.K.q,0)
    error('Only LPs or convex QPs allowed as inner problem');
end

if options.bilevel.rootcuts & (~(isequal(outer_p.K.s,0) & isequal(outer_p.K.q,0)))
    disp('Disjunctive cuts currently only supported when inner is a QP')
    options.bilevel.rootcuts = 0;
end

FRP0 = inner_p;

[merged_mt,merged_vt] = mergemonoms(inner_p,outer_p);

if ~isequal(inner_p.used_variables,outer_p.used_variables)
    invar  = inner_p.used_variables;
    outvar = outer_p.used_variables;
    
    binary_variables = unique([inner_p.used_variables(inner_p.binary_variables) outer_p.used_variables(outer_p.binary_variables)]);
    integer_variables = unique([inner_p.used_variables(inner_p.integer_variables) outer_p.used_variables(outer_p.integer_variables)]);
    semi_variables = unique([inner_p.used_variables(inner_p.semicont_variables) outer_p.used_variables(outer_p.semicont_variables)]);
    all_variables = unique([inner_p.used_variables outer_p.used_variables]);
    if ~isequal(all_variables,inner_p.used_variables )
        inner_p = pad(inner_p,all_variables);
        FRP0 = inner_p;
        FRP0.monomtable = speye(length(inner_p.c));
    end
    if ~isequal(all_variables,outer_p.used_variables )
        outer_p = pad(outer_p,all_variables);
    end
else
    binary_variables = unique([inner_p.used_variables(inner_p.binary_variables) outer_p.used_variables(outer_p.binary_variables)]);
    integer_variables = unique([inner_p.used_variables(inner_p.integer_variables) outer_p.used_variables(outer_p.integer_variables)]);
    semi_variables = unique([inner_p.used_variables(inner_p.semicont_variables) outer_p.used_variables(outer_p.semicont_variables)]);
    all_variables = inner_p.used_variables;
end

outer_p.monomtable   = merged_mt;
outer_p.variabletype = merged_vt;
inner_p.variabletype = merged_vt;
inner_p.monomtable   = merged_vt;
% Index to inner variables
for i = 1:length(y)
    y_var(i) = find(all_variables == getvariables(y(i)));
end

% Index to outer variables
x_var = setdiff(1:length(all_variables),y_var);

% Index to binary variables
bin_var = [];
for i = 1:length(binary_variables)
    bin_var(i) = find(all_variables == binary_variables(i));
end
int_var = [];
for i = 1:length(integer_variables)
    int_var(i) = find(all_variables == integer_variables(i));
end
semi_var = [];
for i = 1:length(semi_variables)
    semi_var(i) = find(all_variables == semi_variables(i));
end

if ~isempty(intersect(y_var,bin_var))
   error('Only LPs or convex QPs allowed as inner problem (inner variables can not be binary)');
end
if ~isempty(intersect(y_var,int_var))
   error('Only LPs or convex QPs allowed as inner problem (inner variables can not be integer)');
end
if ~isempty(intersect(y_var,semi_var))
   error('Only LPs or convex QPs allowed as inner problem (inner variables can not be semi-continuous)');
end
inner_p.binary_variables = bin_var;
outer_p.binary_variables = bin_var;
inner_p.integer_variables = int_var;
outer_p.integer_variables = int_var;
inner_p.semicont_variables = semi_var;
outer_p.semicont_variables = semi_var;

% Number of inequalities in inner model = #bounded dual variables
ninequalities = inner_p.K.l;
nequalities   = inner_p.K.f;

% Add dual related to inequalities in inner model to the model
dual_var = length(all_variables)+1:length(all_variables)+ninequalities;
%dual_var = length(inner_p.c)+1:length(inner_p.c)+ninequalities;

% Add dual related to inequalities in inner model to the model
eqdual_var = dual_var(end)+1:dual_var(end)+inner_p.K.f;

% No cost of duals in outer objective
p = outer_p;
if ~isempty(dual_var)
    p.c(dual_var(end))=0;
    p.Q(dual_var(end),dual_var(end)) = sparse(0);
end

if ~isempty(eqdual_var)
    p.c(eqdual_var(end))=0;
    p.Q(eqdual_var(end),eqdual_var(end)) = sparse(0);
end

% Structure of the constraints
%
%  Stationary equalities
%  Outer equalities
%  Inner equalities
%  Inner LP inequalities
%  Duals positive
%  Outer inequalities (LP, SOCP, SDP)

% Add stationarity to outer model
stationary = [inner_p.c(y_var) 2*inner_p.Q(y_var,:)];
if length(dual_var)>0
    stationary = [stationary -inner_p.F_struc(1+inner_p.K.f:inner_p.K.f+inner_p.K.l,1+y_var)'];
end
if length(eqdual_var)>0
    stationary = [stationary  -inner_p.F_struc(1:inner_p.K.f,1+y_var)'];
end

p.F_struc = [stationary;p.F_struc spalloc(size(p.F_struc,1),length(dual_var) + length(eqdual_var),0)];
p.K.f = p.K.f + length(y_var);

% Add dual>0 to outer model
p.F_struc = [p.F_struc(1:p.K.f,:);spalloc(ninequalities,length(x_var)+length(y_var)+1,0) speye(ninequalities) spalloc(ninequalities,nequalities,0);p.F_struc(1+p.K.f:end,:)];
p.K.l = p.K.l + ninequalities;

% Add inner level constraints to outer model
p.F_struc = [p.F_struc(1:p.K.f,:);inner_p.F_struc spalloc(ninequalities+nequalities,ninequalities+nequalities,0);p.F_struc(1+p.K.f:end,:)];
p.K.f = p.K.f + inner_p.K.f;
p.K.l = p.K.l + inner_p.K.l;
slack_index = p.K.f+1:+p.K.f+ninequalities;

%p.lb = outerinner_p.lb;
%p.ub = outerinner_p.ub;
p.lb(dual_var) = 0;
p.ub(dual_var) =  inf;
p.lb(eqdual_var) = -inf;
p.ub(eqdual_var) =  inf;
p.x0 = [];


%p.variabletype = outerinner_p.variabletype;
%p.monomtable   = outerinner_p.monomtable;
%p.evalMap = outerinner_p.evalMap;
%p.evalVariables = outerinner_p.evalVariables;
for i = 1:length(dual_var)
    p.monomtable(dual_var(i),dual_var(i))=1;
    p.variabletype(dual_var(i)) = 0;
end
for i = 1:length(eqdual_var)
    p.monomtable(eqdual_var(i),eqdual_var(i))=1;
    p.variabletype(eqdual_var(i)) = 0;
end

% xy = sdpvar(length(x_var)+length(y_var),1);
% z = sdpvar(length(dual_var),1);
% res = p.F_struc*[1;xy;z]
% 
% F_bilevel = [res(1:p.K.f) == 0,res(p.K.f+1:end)>0]
% 

% Enable outer problem to be nonconvex etc
p = build_recursive_scheme(p);
% Turned off, generates crash. Unit test in test_bilevel_1
% p = compress_evaluation_scheme(p);

p.lower = -inf;
p.options.verbose = max([0 options.verbose-1]);
p.level = 0;
p.as_free = true(ninequalities,1);
list{1} = p;
lower = -inf;
upper = inf;
iter = 0;
tol = 1e-8;
ndomcuts = 0;
ninfeascuts = 0;

% Extract the inequalities in the inner problem. These are really the
% interesting ones
inner_p.F_struc = [inner_p.F_struc(1+inner_p.K.f:end,:) spalloc(inner_p.K.l,ninequalities+nequalities,0)];

if options.verbose
    disp('* Starting YALMIP bilevel solver.');
    disp(['* Outer solver   : ' outer_p.solver.tag]);
    disp(['* Inner solver   : ' inner_p.solver.tag]);
    disp(['* Max iterations : ' num2str(p.options.bnb.maxiter)]);
    disp(' Node       Upper       Gap(%)      Lower    Open');
end
gap = inf;
xsol = [];
sol.problem = 0;
iter = 0;
inner_p = detectdisjoint(inner_p);

while length(list)>0 & gap > options.bilevel.relgaptol & iter < options.bilevel.maxiter
    iter = iter + 1;
    [p,list,lower] = select(list);

    Comment = '';
    if p.lower<upper

        if strcmp(p.options.solver,'bmibnb') & ~isinf(upper)
            % Allow early termination in bmibnb if it is used in outer
            % probllem
            p.options.bmibnb.lowertarget = upper;
        end
        
        output  = feval(p.solver.call,p);
       
        if output.problem==2
            Comment = 'Unbounded node';
        end
        if output.problem==1
            % Infeasible
            ninfeascuts = ninfeascuts + 1;
            Comment = 'Infeasible in solver';
        else
            
            switch output.problem
                case 0
                    Comment = 'Solved to optimality';
                otherwise
                    Comment = yalmiperror(output.problem);                        
            end

            z = apply_recursive_evaluation(p,output.Primal);
            cost = z'*p.Q*z + p.c'*z + p.f;
                  
            ActuallyFeasible = checkfeasiblefast(p,z,options.bilevel.feastol);
            if ~ActuallyFeasible               
                % Try to solve the relaxed feasibility problem using the
                % inner solver (i.e. treat as LP. If infeasible, it is for
                % sure infeasible
                pAux = p;pAux.c = p.c*0;pAux.Q = p.Q*0;
                outputCheck = feval(inner_p.solver.call,pAux);
                if outputCheck.problem == 1
                    % We will not continue branching, and let the user now
                    % that this choice
                    Comment = ['Infeasible'];
                    cost = inf;
                    sol.problem = 4;  
                else
                    % Hard to say anything                    
                    Comment = ['Infeasible solution returned, resolve => continue'];
                    sol.problem = 4;  
                    cost = p.lower;
                end
            end
            if cost<inf
                 
                 if strcmp(p.options.solver,'bmibnb') 
                     if output.problem == -6
                         sol.problem = -6;
                         sol.info = yalmiperror(-6);
                         info = [];
                         return
                     end
                     p.lb = max([p.lb output.extra.propagatedlb],[],2);
                     p.ub = min([p.ub output.extra.propagatedub],[],2);
                 end
                  
                % These are duals in the original inner problem
                lambda = output.Primal(dual_var);
                % Constraint slacks in original inner problem
                slack = inner_p.F_struc*[1;output.Primal];
                % Outer variables
                xi = z(x_var);
                % Inner variables
                yi = z(y_var);

                res = (slack).*lambda;
                if ActuallyFeasible
                    res = (slack).*lambda;
                else
                    % Generate a dummy residual, to make sure we branch on
                    % the first free 
                    res = (slack).*lambda*0;
                    res(find(p.as_free)) = 1:length(find(p.as_free));
                end
                if (all(p.as_free==0) | max(abs(res(p.as_free)))<options.bilevel.compslacktol) & ActuallyFeasible
                    % Feasible!
                    if upper>cost
                        upper = cost;
                        xsol = xi;
                        zsol = yi;
                        dualsol = output.Primal(dual_var);
                    end
                elseif cost>upper-1e-10
                    ndomcuts = ndomcuts + 1;
                else

                    % No official code, just playing around
                    if  ActuallyFeasible & options.bilevel.solvefrp
                        FRP = FRP0;
                        if 0
                            FRP = fixvariables(FRP0,x_var,xi,y_var);
                        else
                            FRP.F_struc = [xi -sparse(1:length(x_var),x_var,ones(length(x_var),1),length(x_var),length(x_var)+length(y_var));FRP.F_struc];
                            FRP.K.f = FRP.K.f + length(xi);
                            FRP.options.verbose = 0;
                            QQ = sparse(FRP0.Q);
                            cc = sparse(FRP0.c);
                            FRP.c(y_var) = FRP.c(y_var) + 2*FRP.Q(x_var,y_var)'*xi;
                            FRP.Q(x_var,y_var)=0;
                            FRP.Q(y_var,x_var)=0;
                            FRP.Q(x_var,x_var)=0;
                        end
                        outputFRP = feval(inner_p.solver.call,FRP);
                        if outputFRP.problem == 0
                            if 0
                                z = zeros(length(outer_p.c),1);
                                z(x_var) = xi;
                                z(y_var) = outputFRP.Primal;
                                z2 = apply_recursive_evaluation(p,z);
                            else
                                z2 = apply_recursive_evaluation(p,outputFRP.Primal);
                            end
                            costFRP = z2'*outer_p.Q*z2 + outer_p.c'*z2 + outer_p.f;
                            if costFRP < upper & isfeasible(outer_p,z2)
                                upper = costFRP;
                                xsol = z2(x_var);
                                zsol = z2(y_var);
                            end
                        end
                    end

                    [ii,jj_tmp] = max(res(p.as_free));
                    ind_tmp = (1:length(res))';
                    ind_tmp = ind_tmp(p.as_free);
                    jj = ind_tmp(jj_tmp); 
                    
                    if strcmp(p.options.solver,'bmibnb')
                        % Since BMIBNB solves a relaxation of relaxation, it
                        % can generate a lower bound which is lower than
                        % the lower bound before a compl. slack constraint
                        % was added. 
                        p.lower = max(output.lower,lower);
                    else
                        p.lower = cost;
                    end

                    if iter<=options.bilevel.rootcuts
                        % Add a disjunction cut
                        p = disjunction(p,dual_var(jj),inner_p.F_struc(jj,:),output.Primal);
                        % Put in queuee, it will be pulled back immediately
                        list = {list{:},p};
                    else

                        p1 = p;
                        p2 = p;

                        % Add dual == 0 on p1
                        p1.K.f = p1.K.f + 1;
                        p1.F_struc = [zeros(1,size(p1.F_struc,2));p1.F_struc];
                        p1.F_struc(1,1+dual_var(jj))=1;
                        p1.lb(dual_var(jj)) = -inf;
                        p1.ub(dual_var(jj)) = inf;
                        newequality = p1.F_struc(1,:);
                        redundantinequality = findrows(p1.F_struc(p1.K.f+1:end,:),newequality);
                        if ~isempty(redundantinequality)
                            p1.F_struc(p1.K.f+redundantinequality,:)=[];
                            p1.K.l = p1.K.l-length(redundantinequality);
                        end

                        % Add slack == 0
                        p2.K.f = p2.K.f + 1;
                        newequality = inner_p.F_struc(jj,:);
                        p2.F_struc = [newequality;p2.F_struc];
                        redundantinequality = findrows(p2.F_struc(p2.K.f+1:end,:),newequality);
                        if ~isempty(redundantinequality)
                            p2.F_struc(p2.K.f+redundantinequality,:)=[];
                            p2.K.l = p2.K.l-length(redundantinequality);
                        end
                        
                        p1.as_free(jj) = false;
                        p2.as_free(jj) = false;
                        if ~isempty(inner_p.disjoints)
                            here = find(inner_p.disjoints(:,1) == j);
                            if ~isempty(here)
                                p1.as_free(inner_p.disjoints(here,2))=false;
                                p2.as_free(inner_p.disjoints(here,2))=false;
                            else
                                here = find(inner_p.disjoints(:,2) == j);
                                if ~isempty(here)
                                    p1.as_free(inner_p.disjoints(here,1))=false;
                                    p2.as_free(inner_p.disjoints(here,1))=false;
                                end
                            end
                        end
                        
                        p1.level = p.level+1;
                        p2.level = p.level+1;
                        list = {list{:},p1};
                        list = {list{:},p2};
                    end
                end
            end
        end
    else
        ndomcuts = ndomcuts + 1;
    end

    [list,lower] = prune(list,upper);

    gap = abs((upper-lower)/(1e-3+abs(upper)+abs(lower)));
    if isnan(gap)
        gap = inf;
    end
    if options.verbose
        fprintf(' %4.0f : %12.3E  %7.2f   %12.3E  %2.0f  %s\n',iter,full(upper),100*full(gap),full(lower),length(list),Comment)
    end
end
info.upper = upper;
info.iter = iter;
info.ninfeascuts = ninfeascuts;
info.ndomcuts = ndomcuts;

if ~isempty(xsol)
    assign(recover(all_variables(x_var)),xsol);
    assign(recover(all_variables(y_var)),zsol);
else
    sol.problem = 1;
end




function [list,lower] = prune(list,upper)

l = [];
for i = 1:length(list)
    l = [l list{i}.lower];
end
j = find(upper > l+1e-10);
list = {list{j}};
if length(list) == 0
    lower = upper;
else
    lower = min(l(j));
end




function [p,list,lower] = select(list)
l = [];
for i = 1:length(list)
    l = [l list{i}.lower];
end
[i,j] = min(l);

p = list{j};
list = {list{1:j-1},list{j+1:end}};
lower = min(l);

function p = addzero(p,i);

p.K.f = p.K.f + 1;
p.F_struc = [zeros(1,size(p.F_struc,2));p.F_struc];
p.F_struc(1,1+i)=1;


function outer_p = pad(outer_p,all_variables)

[i,loc] = find(ismember(all_variables,outer_p.used_variables));
p = outer_p;
% Set all bounds to infinite, and then place the known bounds
p.lb = -inf(length(all_variables),1);
p.lb(loc) = outer_p.lb;
p.ub = inf(length(all_variables),1);
p.ub(loc) = outer_p.ub;

% Set all variables as linear
p.variabletype = zeros(1,length(all_variables));
p.variabletype(loc) = outer_p.variabletype;

p.c = spalloc(length(all_variables),1,0);
p.c(loc) = outer_p.c;

if ~isempty(p.F_struc)
    p.F_struc = spalloc(size(p.F_struc,1),length(all_variables)+1,nnz(p.F_struc));
    p.F_struc(:,1) = outer_p.F_struc(:,1);
    p.F_struc(:,1+loc) = outer_p.F_struc(:,2:end);
end

% if ~isempty(p.binary_variables)
% end

p.Q = spalloc(length(all_variables),length(all_variables),nnz(outer_p.Q));
p.Q(loc,loc) = outer_p.Q;

outer_p = p;



function p = disjunction(p,variable,const,xstar)
neq = p.K.f+1;

x = sdpvar(length(p.c),1);
e = p.F_struc*[1;x];
Model1 = [x(variable)==0,-e(1:p.K.f)==0, e(1+p.K.f:end)>=0];
Model2 = [const*[1;x]==0,-e(1:p.K.f)==0, e(1+p.K.f:end)>=0];
Ab1 = getbase(sdpvar(Model1));
Ab2 = getbase(sdpvar(Model2));
b1 = -Ab1(:,1);
A1 =  Ab1(:,2:end);
b2 = -Ab2(:,1);
A2 = Ab2(:,2:end);

% b1c = [0;-p.F_struc(:,1)];
% b2c = [const(1);-p.F_struc(:,1)];
% A1c = [-eyev(length(p.c),variable)';p.F_struc(:,2:end)];
% A2c = [-const(2:end);p.F_struc(:,2:end)];
%norm(b1-b1c)
%norm(b2-b2c)
%norm(A1-A1c,inf)
%norm(A2-A2c,inf)

alpha = sdpvar(length(xstar),1);
beta = sdpvar(1);
mu1 = sdpvar(length(b1),1);
mu2 = sdpvar(length(b2),1);

Objective  = alpha'*xstar-beta;
Constraint = [alpha' == mu1'*A1,alpha' == mu2'*A2,beta <= mu1'*b1, beta <= mu2'*b2,mu1(neq+1:end)>0,mu2(neq+1:end)>0];
%Constraint = [alpha' == mu1'*A1,alpha' == mu2'*A2,beta == mu1'*b1, beta == mu2'*b2,mu1(neq+1:end)>0,mu2(neq+1:end)>0];
%Constraint = [Constraint,-10<alpha<10,sum(mu1(neq+1:end))-sum(mu1(1:neq))<10,sum(mu2(neq+1:end))-sum(mu2(1:neq))<10];
%Constraint = [Constraint,-1<alpha<1,mu1(1)+mu2(1) == 1];
Constraint = [Constraint,-1<alpha<1,sum(mu1)+sum(mu2)==1];
%Constraint = [Constraint,sum(mu1(neq+1:end))-sum(mu1(1:neq))<10,sum(mu2(neq+1:end))-sum(mu2(1:neq))<10];

solvesdp(Constraint,Objective,sdpsettings('verbose',0));
p.K.l = p.K.l + 1;
p.F_struc = [p.F_struc;-double(beta) double(alpha)'];


function p = disjunctionFAST(p,variable,const,xstar)
neq = p.K.f+1;
n = length(p.c);

b1 = [0;-p.F_struc(:,1)];
b2 = [const(1);-p.F_struc(:,1)];
A1 = [-eyev(length(p.c),variable)';p.F_struc(:,2:end)];
A2 = [-const(2:end);p.F_struc(:,2:end)];

alpha_ind = 1:length(xstar);
beta_ind = alpha_ind(end)+1;
mu1_ind = (1:length(b1))+beta_ind;
mu2_ind = (1:length(b2))+mu1_ind(end);

alpha = sdpvar(length(xstar),1);
beta = sdpvar(1);
mu1 = sdpvar(length(b1),1);
mu2 = sdpvar(length(b2),1);

p_hull = p;
p_hull.c = zeros(mu2_ind(end),1);
p_hull.c(alpha_ind) = xstar;
p_hull.c(beta_ind) = -1;

% equalities alpha = Ai'*mui, sum(mu)==1
p_hull.K.f = length(xstar)*2+1;
p_hull.F_struc = [zeros(length(xstar),1) eye(length(xstar)) zeros(length(xstar),1) -A1' zeros(length(xstar),length(b2))];
p_hull.F_struc = [p_hull.F_struc;
    zeros(length(xstar),1) eye(length(xstar)) zeros(length(xstar),1) zeros(length(xstar),length(b2)) -A2'];
p_hull.F_struc = [p_hull.F_struc ;1 zeros(1,length(xstar)) 0 -ones(1,length(b1)+length(b2))];

% Inequalities
p_hull.F_struc = [p_hull.F_struc ;0 zeros(1,length(xstar)) -1 b1' b2'*0];
p_hull.F_struc = [p_hull.F_struc ;0 zeros(1,length(xstar)) -1 0*b1' b2'];
npmu = length(b1)-neq;
p_hull.F_struc = [p_hull.F_struc; zeros(npmu,1) zeros(npmu,length(xstar)) zeros(npmu,1) zeros(npmu,neq) eye(npmu) zeros(npmu,length(b2))];
p_hull.F_struc = [p_hull.F_struc; zeros(npmu,1) zeros(npmu,length(xstar)) zeros(npmu,1) zeros(npmu,length(b1)) zeros(npmu,neq) eye(npmu)];
p_hull.F_struc = [p_hull.F_struc; ones(length(xstar),1) -eye(length(xstar)) zeros(length(xstar),1+2*length(b1))];
p_hull.F_struc = [p_hull.F_struc; ones(length(xstar),1)  eye(length(xstar)) zeros(length(xstar),1+2*length(b1))];
p_hull.K.l = 2+npmu*2+2*length(xstar);
p_hull.lb = [];
p_hull.ub = [];
output = feval(p_hull.solver.call,p_hull);
alpha = output.Primal(alpha_ind);
beta = output.Primal(beta_ind);

% Objective  = alpha'*xstar-beta;
% Constraint = [alpha' == mu1'*A1,alpha' == mu2'*A2,beta <= mu1'*b1, beta <= mu2'*b2,mu1(neq+1:end)>0,mu2(neq+1:end)>0];
% Constraint = [Constraint,-1<alpha<1,sum(mu1)+sum(mu2)==1];
% Constraint = [Constraint,sum(mu1)+sum(mu2)==1];
% solvesdp(Constraint,Objective,sdpsettings('verbose',0));

p.K.l = p.K.l + 1;
p.F_struc = [p.F_struc;-double(beta) double(alpha)'];


function feas = isfeasible(p,x)

feas = checkfeasiblefast(p,x,1e-8);


function p = detectdisjoint(p);

p.disjoints = [];
% for i = 1:p.K.l
%     row1 = p.F_struc(i+p.K.f,:);
%     for j = 2:1:p.K.l
%         row2 = p.F_struc(j+p.K.f,:);
%         
%         if all(abs(row1)-abs(row2)==0)
%             % candidate
%             if nnz(row1 == -row2 & row1~=0)==1
%                 p.disjoints = [p.disjoints;i j];
%             end
%         end
%     end
% end
% 


function FRP = fixvariables(FRP0,x_var,xi,y_var);
% Copy current model                        
FRP = FRP0;

% 
FRP.c(y_var) = FRP.c(y_var) + 2*FRP.Q(x_var,y_var)'*xi;
FRP.c(x_var) = [];
FRP.Q(:,x_var) = [];
FRP.Q(x_var,:) = [];
FRP.lb(x_var) = [];
FRP.ub(x_var) = [];
B = FRP.F_struc(:,1+x_var);
FRP.F_struc(:,1+x_var)=[];
FRP.F_struc(:,1) = FRP.F_struc(:,1) + B*xi;

%                         FRP.F_struc = [xi -sparse(1:length(x_var),x_var,ones(length(x_var),1),length(x_var),length(x_var)+length(y_var));FRP.F_struc];
%                         FRP.K.f = FRP.K.f + length(xi);
%                         FRP.options.verbose = 0;
%                         QQ = FRP0.Q;
%                         cc = FRP0.c;
%                         FRP.c(y_var) = FRP.c(y_var) + 2*FRP.Q(x_var,y_var)'*xi;
%                         FRP.Q(x_var,y_var)=0;
%                         FRP.Q(y_var,x_var)=0;
%                         FRP.Q(x_var,x_var)=0;



function [merged_mt,merged_vt] = mergemonoms(inner_p,outer_p);

if isequal(inner_p.used_variables,outer_p.used_variables)
    merged_mt = inner_p.monomtable;
    merged_vt = inner_p.variabletype;
else
    invar  = inner_p.used_variables;
    outvar = outer_p.used_variables;
    all_variables = unique([invar outvar]);
    
    [i_inner,loc_inner] = find(ismember(all_variables,inner_p.used_variables));
    [i_outer,loc_outer] = find(ismember(all_variables,outer_p.used_variables));

    merged_mt = spalloc(length(all_variables),length(all_variables),0);
    merged_vt = zeros(1,length(all_variables));
    
    for i = 1:length(i_inner)
        [ii,jj,kk] = find(inner_p.monomtable(i,:));
        merged_mt(loc_inner(i),loc_inner(jj)) = kk;
        merged_vt(loc_inner(i)) = inner_p.variabletype(i);
    end
    for i = 1:length(i_outer)
        [ii,jj,kk] = find(outer_p.monomtable(i,:));
        merged_mt(loc_outer(i),loc_outer(jj)) = kk;
        merged_vt(loc_outer(i)) = outer_p.variabletype(i);
    end   
end
    
    
    

