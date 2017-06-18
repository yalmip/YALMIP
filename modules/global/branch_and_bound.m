function [x_min,solved_nodes,lower,upper,lower_hist,upper_hist,timing,counter] = branch_and_bound(p,x_min,upper,timing)

% *************************************************************************
% Create handles to solvers
% *************************************************************************
lowersolver = p.solver.lowersolver.call; % For relaxed lower bound problem
uppersolver = p.solver.uppersolver.call; % Local nonlinear upper bound
lpsolver    = p.solver.lpsolver.call;    % LP solver for bound propagation

% *************************************************************************f
% GLOBAL PROBLEM DATA (these variables are the same in all nodes)
% *************************************************************************
c       = p.c;
Q       = p.Q;
f       = p.f;
K       = p.K;
options = p.options;

% *************************************************************************
% DEFINE UPPER BOUND PROBLEM. Basically just remove the cuts
% *************************************************************************
p_upper = cleanuppermodel(p);

% *************************************************************************
% Active constraints in main model
% 0   : Inactive constraint (i.e. a cut which unused)
% 1   : Active constraint
% inf : Removed constraint  (found to be redundant)
% *************************************************************************
p.InequalityConstraintState = ones(p.K.l,1);
p.InequalityConstraintState(p.KCut.l,1) = 0;
p.EqualityConstraintState = ones(p.K.f,1);

% *************************************************************************
% LPs ARE USED IN  BOX-REDUCTION
% *************************************************************************
p.lpcuts = p.F_struc(1+p.K.f:1:p.K.l+p.K.f,:);
p.cutState = ones(p.K.l,1);
p.cutState(p.KCut.l,1) = 0; % Don't use to begin with

% *************************************************************************
% INITIALITAZION
% *************************************************************************
p.depth = 0;        % depth in search tree
p.dpos  = 0;        % used for debugging
p.lower = NaN;
lower   = NaN;
gap     = inf;
stack   = [];
solved_nodes = 0;
numGlobalSolutions = 0;

% *************************************************************************
% Silly hack to speed up solver calls
% *************************************************************************
p.getsolvertime = 0;

counter = p.counter;

if options.bmibnb.verbose>0
    disp('* Starting YALMIP global branch & bound.');
    disp(['* Branch-variables : ' num2str(length(p.branch_variables))]);
    disp(['* Upper solver     : ' p.solver.uppersolver.tag]);
    disp(['* Lower solver     : ' p.solver.lowersolver.tag]);
    if p.options.bmibnb.lpreduce
        disp(['* LP solver        : ' p.solver.lpsolver.tag]);
    end
    disp(' Node       Upper      Gap(%)       Lower    Open');
end

t_start = cputime;
go_on  = 1;

reduction_result = [];
lower_hist = [];
upper_hist = [];
p.branchwidth = [];

pseudo_costgain=[];
pseudo_variable=[];

while go_on

    % *********************************************************************
    % ASSUME THAT WE WON'T FATHOME
    % *********************************************************************
    keep_digging = 1;

    % *********************************************************************
    % Strenghten variable bounds a couple of runs
    % *********************************************************************
    p.changedbounds = 1;
    
    for i = 1:length(options.bmibnb.strengthscheme)
        if ~p.feasible
            break
        end      
        switch options.bmibnb.strengthscheme(i)
            case 1
                p = updatebounds_recursive_evaluation(p);
            case 2
                p = updateboundsfromupper(p,upper,p.originalModel);
            case 3
                p = propagatequadratics(p);
            case 4
                p = propagate_bounds_from_complementary(p);
            case 5
                tstart = tic;
                p = domain_reduction(p,upper,lower,lpsolver,x_min);
                timing.domainreduce = timing.domainreduce + toc(tstart);
            case 6
                p = propagate_bounds_from_equalities(p);
            otherwise
        end
    end

    % *********************************************************************
    % Detect redundant constraints
    % *********************************************************************
    p = remove_redundant(p);

    % *********************************************************************
    % SOLVE LOWER AND UPPER
    % *********************************************************************
    if p.feasible

        [output,cost,p,timing] = solvelower(p,options,lowersolver,x_min,upper,timing);

        if output.problem == -1
            % We have no idea what happened. 
            % Behave as if it worked, so we can branch as see if things
            % clean up nicely
            cost = p.lower;
            if isnan(cost)
                cost = -inf;
            end
            output.problem = 3;
        end
        
        % Cplex sucks...
        if output.problem == 12
            pp = p;
            pp.c = pp.c*0;
            [output2,cost2] = solvelower(pp,options,lowersolver,[],[],timing);
            if output2.problem == 0
                output.problem = 2;
            else
                output.problem = 1;
            end
        end
        
        % GLPK sucks in st_e06
        if abs(p.lb(p.linears)-p.ub(p.linears)) <= 1e-3 & output.problem==1
            x = (p.lb+p.ub)/2;
            z = evaluate_nonlinear(p,x);
            oldCount = numGlobalSolutions;
            if numGlobalSolutions < p.options.bmibnb.numglobal
                [upper,x_min,cost,info_text,numGlobalSolutions] = heuristics_from_relaxed(p_upper,x,upper,x_min,cost,numGlobalSolutions);               
            end
        end

        info_text = '';
        switch output.problem
            case {1,12} % Infeasible
                info_text = 'Infeasible';
                keep_digging = 0;
                cost = inf;
                feasible = 0;

            case 2 % Unbounded (should not happen!)
                cost = -inf;
                x = output.Primal;

            case {0,3,4} % No problems (disregard numerical problems)

                if (output.problem == 3) | (output.problem == 4)
                    info_text = 'Numerical problems in lower bound solver';
                end
                x = output.Primal;

                if ~isempty(p.branchwidth)
                    if ~isempty(p.lower)
                        pseudo_costgain = [pseudo_costgain (cost-p.lower)/p.branchwidth];
                        pseudo_variable = [pseudo_variable p.spliton];
                    end
                end
                % UPDATE THE LOWER BOUND
                if isnan(lower)
                    lower = cost;
                end
                if ~isempty(stack)
                    lower = min(cost,min([stack.lower]));
                else
                    lower = min(lower,cost);
                end

                relgap = 100*(upper-lower)/(1+abs(upper));
                relgap_too_big = (isinf(lower) | isnan(relgap) | relgap>options.bmibnb.relgaptol);
                if cost<upper-1e-5 & relgap_too_big

                    z = evaluate_nonlinear(p,x);

                    % Manage cuts etc
                    p = addsdpcut(p,z);
                    p = addlpcuts(p,x);

                    oldCount = numGlobalSolutions;
                    if numGlobalSolutions < p.options.bmibnb.numglobal                        
                        [upper,x_min,cost,info_text2,numGlobalSolutions] = heuristics_from_relaxed(p_upper,x,upper,x_min,cost,numGlobalSolutions);
                        if length(info_text)==0
                            info_text = info_text2;
                        elseif  length(info_text2)>0
                            info_text = [info_text ' | ' info_text2];
                        else
                            info_text = info_text; 
                        end
                        if ~isequal(p.solver.uppersolver.tag,'none')
                            if upper > p.options.bmibnb.target
                                if options.bmibnb.lowertarget > lower                                    
                                    [upper,x_min,info_text,numGlobalSolutions,timing] = solve_upper_in_node(p,p_upper,x,upper,x_min,uppersolver,info_text,numGlobalSolutions,timing);
                                    p.counter.uppersolved = p.counter.uppersolved + 1;
                                end
                            end
                        end
                    end
                else
                    keep_digging = 0;
                    info_text = 'Poor bound';
                end
            otherwise
                cost = -inf;
                x = (p.lb+p.ub)/2;
        end
    else
        info_text = 'Infeasible';
        keep_digging = 0;
        cost = inf;
        feasible = 0;
    end
    solved_nodes = solved_nodes+1;

    % ************************************************
    % PRUNE SUBOPTIMAL REGIONS BASED ON UPPER BOUND
    % ************************************************
    if ~isempty(stack)
        [stack,lower] = prune(stack,upper,options,solved_nodes,p);
    end
    if isempty(stack)
        if isinf(cost)
            lower = upper;
        else
            lower = cost;
        end
    else
        lower = min(lower,cost);
    end

    % ************************************************
    % CONTINUE SPLITTING?
    % ************************************************
    if keep_digging & max(p.ub(p.branch_variables)-p.lb(p.branch_variables))>options.bmibnb.vartol
        node = [];
      %  already_tested = []
      %  while ~isempty(setdiff(p.branch_variables,already_tested)) & isempty(node)
      %  temp = p.branch_variables;
      %  p.branch_variables=setdiff(p.branch_variables,already_tested);
        spliton = branchvariable(p,options,x);
      %  p.branch_variables = union(p.branch_variables,already_tested);
      %  already_tested = [already_tested spliton];
      
        if ismember(spliton,p.complementary)
            i = find(p.complementary(:,1) == spliton);
            if isempty(i)
                i = find(p.complementary(:,2) == spliton);
            end
            % Either v1 or v2 is zero
            v1 = p.complementary(i,1);
            v2 = p.complementary(i,2);
            gap_over_v1 = (p.lb(v1)<=0) & (p.ub(v1)>=0) & (p.ub(v1)-p.lb(v2))>0;
            gap_over_v2 = (p.lb(v2)<=0) & (p.ub(v2)>=0) & (p.ub(v2)-p.lb(v2))>0;
            
            if gap_over_v1
                pp = p;
                pp.complementary( find((pp.complementary(:,1)==v1) | (pp.complementary(:,2)==v1)),:)=[];
                node = savetonode(pp,v1,0,0,-1,x,cost,p.EqualityConstraintState,p.InequalityConstraintState,p.cutState);
                node.bilinears = p.bilinears;
                node = updateonenonlinearbound(node,spliton);
                if all(node.lb <= node.ub)
                    node.branchwidth=[];
                    stack = push(stack,node);                 
                end
            end
            if gap_over_v2
                pp = p;
                %pp.complementary(i,:)=[];
                pp.complementary( find((pp.complementary(:,1)==v2) | (pp.complementary(:,2)==v2)),:)=[];
                node = savetonode(pp,v2,0,0,-1,x,cost,p.EqualityConstraintState,p.InequalityConstraintState,p.cutState);
                node.bilinears = p.bilinears;
                node = updateonenonlinearbound(node,spliton);
                if all(node.lb <= node.ub)
                    node.branchwidth=[];
                    stack = push(stack,node);                 
                end
            end     
        end
        if isempty(node)
            bounds  = partition(p,options,spliton,x);
            if length(bounds)>3
                error('REPORT BOUND LENGTH UNIMPLEMENTED BUG')
            end
            for i = 1:length(bounds)-1
                if ismember(spliton,union(p.binary_variables,p.integer_variables)) & (i==2)
                    node = savetonode(p,spliton,bounds(i)+1,bounds(i+1),-1,x,cost,p.EqualityConstraintState,p.InequalityConstraintState,p.cutState);
                else
                    node = savetonode(p,spliton,bounds(i),bounds(i+1),-1,x,cost,p.EqualityConstraintState,p.InequalityConstraintState,p.cutState);
                end
                node.bilinears = p.bilinears;
                node = updateonenonlinearbound(node,spliton);
                node.branchwidth = [p.ub(spliton)-p.lb(spliton)];
                if all(node.lb <= node.ub)
                    stack = push(stack,node);
                end
            end
        end
        lower = min([stack.lower]);
    end

    if ~isempty(p)
        counter = p.counter;
    end
    % ************************************************
    %  Pick and create a suitable node
    % ************************************************    
    [p,stack] = selectbranch(p,options,stack,x_min,upper);
    
    if isempty(p)
        if ~isinf(upper)
            relgap = 0;
        end
        if isinf(upper) & isinf(lower)
            relgap = inf;
        end
        depth = 0;
    else
        relgap = 100*(upper-lower)/(1+max(abs(lower)+abs(upper))/2);
        depth = p.depth;
    end
    if options.bmibnb.verbose>0
        fprintf(' %4.0f : %12.3E  %7.2f   %12.3E  %2.0f  %s  \n',solved_nodes,upper,relgap,lower,length(stack)+length(p),info_text);
    end

    absgap = upper-lower;
    % ************************************************
    % Continue?
    % ************************************************
    time_ok = cputime-t_start < options.bmibnb.maxtime;
    iter_ok = solved_nodes < options.bmibnb.maxiter;
    any_nodes = ~isempty(p);
    relgap_too_big = (isinf(lower) | isnan(relgap) | relgap>100*options.bmibnb.relgaptol);
    absgap_too_big = (isinf(lower) | isnan(absgap) | absgap>options.bmibnb.absgaptol);
    uppertarget_not_met = upper > options.bmibnb.target;
    lowertarget_not_met = lower < options.bmibnb.lowertarget;
    go_on = uppertarget_not_met & lowertarget_not_met & time_ok & any_nodes & iter_ok & relgap_too_big & absgap_too_big;
    lower_hist = [lower_hist lower];
    upper_hist = [upper_hist upper];
end
if options.bmibnb.verbose>0
    fprintf(['* Finished.  Cost: ' num2str(upper) ' Gap: ' num2str(relgap) '\n']);
end

%save dummy x_min

% *************************************************************************
% Stack functionality
% *************************************************************************
function stack = push(stackin,p)
if ~isempty(stackin)
    stack = [p;stackin];
else
    stack(1)=p;
end

function [p,stack] = pull(stack,method,x_min,upper,branch_variables);
if ~isempty(stack)
    switch method
        case 'maxvol'
            for i = 1:length(stack)
                vol(i) = sum(stack(i).ub(branch_variables)-stack(i).lb(branch_variables));
            end
            [i,j] = max(vol);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);

        case 'best'
            [i,j]=min([stack.lower]);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);

        otherwise
    end
else
    p =[];
end

function [stack,lower] = prune(stack,upper,options,solved_nodes,p)
if ~isempty(stack)
    toolarge = find([stack.lower]>upper*(1+1e-4));
    if ~isempty(toolarge)
        stack(toolarge)=[];
    end
    if ~isempty(stack)
        
        for j = 1:length(stack)
            if nnz(p.c.*(stack(j).ub-stack(j).lb)) == 1 & nnz(p.Q)==0
                i = find(p.c.*(stack(j).ub-stack(j).lb));
                if p.c(i)>0
                    stack(j).ub(i) = min([stack(j).ub(i) upper]);
                end
            end
        end
        
        indPOS = find(p.c>0);
        indNEG = find(p.c<0);
        LB = [stack.lb];
        UB = [stack.ub];
        LOWER =  p.c([indPOS(:);indNEG(:)])'*[LB(indPOS,:);UB(indNEG,:)];
        toolarge = find(LOWER > upper*(1-1e-8));
        stack(toolarge)=[];
    end
end
if ~isempty(stack)
    lower = min([stack.lower]);
else
    lower = upper;
end

function node = savetonode(p,spliton,bounds1,bounds2,direction,x,cost,EqualityConstraintState,InequalityConstraintState,cutState);
node.lb = p.lb;
node.ub = p.ub;
node.lb(spliton) = bounds1;
node.ub(spliton) = bounds2;
node.lb(p.integer_variables) = ceil(node.lb(p.integer_variables));
node.ub(p.integer_variables) = floor(node.ub(p.integer_variables));
node.lb(p.binary_variables) = ceil(node.lb(p.binary_variables));
node.ub(p.binary_variables) = floor(node.ub(p.binary_variables));
node.complementary = p.complementary;

if direction == -1
    node.dpos = p.dpos-1/(2^sqrt(p.depth));
else
    node.dpos = p.dpos+1/(2^sqrt(p.depth));
end
node.spliton = spliton;
node.depth = p.depth+1;
node.x0 = x;
node.lpcuts = p.lpcuts;
node.lower = cost;
node.InequalityConstraintState = InequalityConstraintState;
node.EqualityConstraintState = EqualityConstraintState;
node.cutState = cutState;

% *************************************
% DERIVE LINEAR CUTS FROM SDPs
% *************************************
function p = addsdpcut(p,x)
if p.K.s > 0
    top = p.K.f+p.K.l+1;
    newcuts = 1;
    newF = [];
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        X = p.F_struc(top:top+n^2-1,:)*[1;x];
        X = full(reshape(X,n,n));
        [d,v] = eig(X);
        for m = 1:length(v)
            if v(m,m)<0
                for j = 1:length(x)+1;
                    newF(newcuts,j)= d(:,m)'*reshape(p.F_struc(top:top+n^2-1,j),n,n)*d(:,m);
                end
                % max(abs(newF(:,2:end)),[],2)
                newF(newcuts,1)=newF(newcuts,1)+1e-6;
                newcuts = newcuts + 1;
                if size(p.lpcuts,1)>0
                    dist = p.lpcuts*newF(newcuts-1,:)'/(newF(newcuts-1,:)*newF(newcuts-1,:)');
                    if any(abs(dist-1)<1e-3)
                        newF = newF(1:end-1,:);
                        newcuts = newcuts - 1;
                    end
                end
            end
        end
        top = top+n^2;
    end

    if ~isempty(newF)
        % Don't keep all
        m = size(newF,2);
        %  size(p.lpcuts)
        p.lpcuts = [newF;p.lpcuts];
        p.cutState = [ones(size(newF,1),1);p.cutState];
        violations = p.lpcuts*[1;x];
        p.lpcuts = p.lpcuts(violations<0.1,:);

        if size(p.lpcuts,1)>15*m
            disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            violations = p.lpcuts*[1;x];
            [i,j] = sort(violations);
            %p.lpcuts = p.lpcuts(j(1:15*m),:);
            %p.cutState = lpcuts = p.lpcuts(j(1:15*m),:);
            %p.lpcuts = p.lpcuts(end-15*m+1:end,:);
        end
    end
end

function p = addlpcuts(p,z)
if ~isempty(p.lpcuts)
    inactiveCuts = find(~p.cutState);
    violation = p.lpcuts(inactiveCuts,:)*[1;z];
    need_to_add = find(violation < -1e-4);
    if ~isempty(need_to_add)
        p.cutState(inactiveCuts(need_to_add)) = 1;
    end
    inactiveCuts = find(p.InequalityConstraintState == 0 );
    violation = p.F_struc(p.K.f+inactiveCuts,:)*[1;z];
    need_to_add = find(violation < -1e-4);
    if ~isempty(need_to_add)
        p.InequalityConstraintState(inactiveCuts(need_to_add)) = 1;
    end
end

% *************************************************************************
% Strategy for deciding which variable to branch on
% *************************************************************************
function spliton = branchvariable(p,options,x)
% Split if box is to narrow
width = abs(p.ub(p.branch_variables)-p.lb(p.branch_variables));
% if ~isempty(p.binary_variables)
%     width_bin = min([abs(1-x(p.binary_variables)) abs(x(p.binary_variables))],[],2);
% end
if isempty(p.bilinears) | ~isempty(p.evalMap) | any(p.variabletype > 2)%(min(width)/max(width) < 0.1) | (size(p.bilinears,1)==0) %
    [i,j] = max(width);%.*p.weight(p.branch_variables));
    spliton = p.branch_variables(j);
else
    res = x(p.bilinears(:,1))-x(p.bilinears(:,2)).*x(p.bilinears(:,3));
    [ii,jj] = sort(abs(res));
    v1 = p.bilinears(jj(end),2);
    v2 = p.bilinears(jj(end),3);

    acc_res1 = sum(abs(res(find((p.bilinears(:,2)==v1) |  p.bilinears(:,3)==v1))));
    acc_res2 = sum(abs(res(find((p.bilinears(:,2)==v2) |  p.bilinears(:,3)==v2))));
    
    if abs(acc_res1-acc_res2)<1e-3 & ismember(v2,p.branch_variables) & ismember(v1,p.branch_variables)
        if abs(p.ub(v1)-p.lb(v1))>abs(p.ub(v2)-p.lb(v2))
              spliton = v1;
        elseif abs(p.ub(v1)-p.lb(v1))<abs(p.ub(v2)-p.lb(v2))
            spliton = v2;
        else
        % Oops, two with the same impact. To avoid that we keep pruning on
        % a variable that doesn't influence the bounds, we flip a coin on
        % which to branch on
        if rand(1)>0.5
            spliton = v1;
        else
            spliton = v2;
        end
        end
    else
        if (~ismember(v2,p.branch_variables) | (acc_res1>acc_res2)) & ismember(v1,p.branch_variables)
            spliton = v1;
        elseif ismember(v2,p.branch_variables)
            spliton = v2;
        else
            [i,j] = max(width);
            spliton = p.branch_variables(j);
        end
    end   
end

% *************************************************************************
% Strategy for diving the search space
% *************************************************************************
function bounds = partition(p,options,spliton,x_min)
x = x_min;
if isinf(p.lb(spliton))
    %bounds = [p.lb(spliton) x_min(spliton) p.ub(spliton)]
    %return
    p.lb(spliton) = -1e6;
end
if isinf(p.ub(spliton))
    %bounds = [p.lb(spliton) x_min(spliton) p.ub(spliton)]
    %return
    p.ub(spliton) = 1e6;
end

switch options.bmibnb.branchrule
    case 'omega'
        if ~isempty(x_min)
            U = p.ub(spliton);
            L = p.lb(spliton);
            x = x(spliton);
            bounds = [p.lb(spliton) 0.5*max(p.lb(spliton),min(x_min(spliton),p.ub(spliton)))+0.5*(p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)];
        else
            bounds = [p.lb(spliton) (p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)];
        end
    case 'bisect'
        bounds = [p.lb(spliton) (p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)];
    otherwise
        bounds = [p.lb(spliton) (p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)];
end
if isnan(bounds(2)) %FIX
    if isinf(p.lb(spliton))
        p.lb(spliton) = -1e6;
    end
    if isinf(p.ub(spliton))
        p.ub(spliton) = 1e6;
    end
    bounds(2) = (p.lb(spliton)+p.ub(spliton))/2;
end

function [p,stack] = selectbranch(p,options,stack,x_min,upper,cost_improvements)
switch options.bmibnb.branchmethod
    case 'maxvol'
        [node,stack] = pull(stack,'maxvol',x_min,upper,p.branch_variables);
    case 'best'
        [node,stack] = pull(stack,'best',x_min,upper);
    case 'best-estimate'
        [node,stack] = pull(stack,'best-estimate',x_min,upper,[],cost_improvements);

    otherwise
        [node,stack] = pull(stack,'best',x_min,upper);
end
% Copy node data to p
if isempty(node)
    p = [];
else
    p.depth = node.depth;
    p.dpos = node.dpos;
    p.spliton = node.spliton;
    p.lb = node.lb;
    p.ub = node.ub;
    p.lower = node.lower;
    p.lpcuts = node.lpcuts;
    p.x0 = node.x0;
    p.InequalityConstraintState = node.InequalityConstraintState;
    p.EqualityConstraintState = node.EqualityConstraintState;
    p.complementary = node.complementary;
    p.cutState = node.cutState;
    p.feasible = 1;
    p.branchwidth = node.branchwidth;
end

% *************************************************************************
% Heuristics from relaxed
% Basically nothing coded yet. Just check feasibility...
% *************************************************************************
function [upper,x_min,cost,info_text,numglobals] = heuristics_from_relaxed(p_upper,x,upper,x_min,cost,numglobals)
%load dummy;U = [x(1) x(2) x(4);0 x(3) x(5);0 0 x(6)];P=U'*U;i = find(triu(ones(length(A))-eye(length(A))));-log(det(U'*U))+trace(A*U'*U)+2*sum(invsathub(P(i),lambda))
x(p_upper.binary_variables) = round(x(p_upper.binary_variables));
x(p_upper.integer_variables) = round(x(p_upper.integer_variables));

z = apply_recursive_evaluation(p_upper,x(1:length(p_upper.c)));
%z = evaluate_nonlinear(p_upper,x);

relaxed_residual = constraint_residuals(p_upper,z);

eq_ok = all(relaxed_residual(1:p_upper.K.f)>=-p_upper.options.bmibnb.eqtol);
iq_ok = all(relaxed_residual(1+p_upper.K.f:end)>=p_upper.options.bmibnb.pdtol);

relaxed_feasible = eq_ok & iq_ok;
info_text = '';
if relaxed_feasible
    this_upper = p_upper.f+p_upper.c'*z+z'*p_upper.Q*z;
    if (this_upper < (1-1e-5)*upper) & (this_upper < upper - 1e-5)
        x_min = x;
        upper = this_upper;
        info_text = 'Improved solution';
        cost = cost-1e-10; % Otherwise we'll fathome!
        numglobals = numglobals + 1;
    end
end

% *************************************************************************
% Detect redundant constraints
% *************************************************************************
function p = remove_redundant(p);

if isempty(p.F_struc)
    return
end

b = p.F_struc(1+p.K.f:p.K.l+p.K.f,1);
A = -p.F_struc(1+p.K.f:p.K.l+p.K.f,2:end);

redundant = find(((A>0).*A*(p.ub-p.lb) - (b-A*p.lb) <-1e-2));

if length(redundant)>1
    p.InequalityConstraintState(redundant) = inf;
end

if p.options.bmibnb.lpreduce
    b = p.lpcuts(:,1);
    A = -p.lpcuts(:,2:end);
    redundant = find(((A>0).*A*(p.ub-p.lb) - (b-A*p.lb) <-1e-2));
    if length(redundant)>1
        p.lpcuts(redundant,:) = [];
        p.cutState(redundant) = [];
    end
end

if p.K.f > 0
    b = p.F_struc(1:p.K.f,1);
    A = -p.F_struc(1:p.K.f,2:end);
    s1 = ((A>0).*A*(p.ub-p.lb) - (b-A*p.lb) <1e-6);
    s2 = ((-A>0).*(-A)*(p.ub-p.lb) - ((-b)-(-A)*p.lb) <1e-6);
    redundant = find(s1 & s2);
    if length(redundant)>1
        p.EqualityConstraintState(redundant) = inf;
    end
end

% *************************************************************************
% Clean the upper bound model
% Remove cut constraints, and
% possibly unused variables not used
% *************************************************************************
function p = cleanuppermodel(p);

% We might have created a bilinear model from an original polynomial model.
% We should use the original model when we solve the upper bound problem.
p_bilinear = p;
p = p.originalModel;

% Remove cuts
p.F_struc(p.K.f+p.KCut.l,:)=[];
p.K.l = p.K.l - length(p.KCut.l);
n_start = length(p.c);

% Quadratic mode, and quadratic aware solver?
bilinear_variables = find(p.variabletype == 1 | p.variabletype == 2);
if ~isempty(bilinear_variables)
    used_in_c = find(p.c);
    quadraticterms = used_in_c(find(ismember(used_in_c,bilinear_variables)));
    if ~isempty(quadraticterms) & p.solver.uppersolver.objective.quadratic.nonconvex
        usedinquadratic = zeros(1,length(p.c));
        for i = 1:length(quadraticterms)
            Qij = p.c(quadraticterms(i));
            power_index = find(p.monomtable(quadraticterms(i),:));
            if length(power_index) == 1
                p.Q(power_index,power_index) = Qij;
            else
                p.Q(power_index(1),power_index(2)) = Qij/2;
                p.Q(power_index(2),power_index(1)) = Qij/2;
            end
            p.c(quadraticterms(i)) = 0;
        end
    end
end

% Remove SDP cuts
if length(p.KCut.s)>0
    starts = p.K.f+p.K.l + [1 1+cumsum((p.K.s).^2)];
    remove_these = [];
    for i = 1:length(p.KCut.s)
        j = p.KCut.s(i);
        remove_these = [remove_these;(starts(j):starts(j+1)-1)'];
    end
    p.F_struc(remove_these,:)=[];
    p.K.s(p.KCut.s) = [];
end
p.lb = p_bilinear.lb(1:length(p.c));
p.ub = p_bilinear.ub(1:length(p.c));
p.bilinears = [];