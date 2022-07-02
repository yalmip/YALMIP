function [x_min,solved_nodes,lower,upper,profile,diagnostics] = bnb_branch_and_bound(p)

%% *******************************
% We don't need this
% *******************************
p.options.savesolveroutput = 0;
p.options.saveduals = 0;
p.options.dimacs = 0;
diagnostics = 0;
bnbsolvertime = clock;

%% *******************************
% Tracking performance etc
% *******************************
profile.local_solver_time = 0;

%% *************************************************************************
% We save this to re-use some stuff in fmincon
% *************************************************************************
p.options.savesolverinput = 1;

%% *******************************
% SET-UP ROOT PROBLEM
% *******************************
p.depth = 0;
p.lower = NaN;
[p,x_min,upper] = initializesolution(p,[],inf);

%%
% See if there is structure which can be exploited to improve primal
% heuristics. If all objective variables only enter SDPs and LPs through
% positive coefficients (and diagonal in SDP), the we can try a bisection
% on these, if the combinatorial stuff is feasible
% FIXME: Exploit
% p = detectMonotoneSDPObjective(p);

%% *******************************
% Global stuff
% *******************************
lower = NaN;
stack = stackCreate;

%% *******************************
% Create function handle to solver
% *******************************
lowersolver = p.solver.lower.call;
uppersolver = p.options.bnb.uppersolver;

p.corig = p.c;

%% *******************************
% INVARIANT PROBLEM DATA
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
    p.integer_variables = [];
    
    for i = 1:length(c)
        p.c = eyev(n,i);
        output = feval(lowersolver,p);
        if output.problem == 0
            p.lb(i) = max(p.lb(i),output.Primal(i)-1e-3);
        end
        p.c = -eyev(n,i);
        output = feval(lowersolver,p);
        if output.problem == 0
            p.ub(i) = min(p.ub(i),output.Primal(i)+1e-3);
        end
        p.lb(saveBinary) = ceil(p.lb(saveBinary)-1e-3);
        p.ub(saveBinary) = floor(p.ub(saveBinary)+1e-3);
    end
    p.binary_variables = saveBinary;
    p.integer_variables = saveInteger;
    
    p.Q = saveQ;
    p.c = savec;
end

%%************************************************
% Some hacks to speed up solver calls
% Only track solver-time if user wants profile
% ************************************************
p.getsolvertime = p.options.bnb.profile;

%%*******************************
%% DISPLAY HEADER
% *******************************
originalDiscrete = [p.integer_variables(:);p.binary_variables(:)];
originalBinary   = p.binary_variables(:);

if nnz(Q-fix(Q))==0 && (nnz(p.c-fix(p.c))==0) && isequal(p.K.m,0)
    can_use_ceil_lower = all(ismember(find(p.c),originalDiscrete));
    can_use_ceil_lower = can_use_ceil_lower && all(ismember(find(any(Q,2)),originalDiscrete));
else
    can_use_ceil_lower = 0;
end

% if nnz(Q)==0 && nnz(c)==1 && ~any(p.K.m)
%     p.simplecost = 1;
% else
%     p.simplecost = 0;
% end

p = detectSOS(p);
p = detect_knapsack(p);
p = detect_atleast(p); 
%p = presolveSOS(p);
%p = smashFixed(p);
p = propagate_bounds_from_qualities(p);

p.options.allowsmashing = 1;
poriginal = p;
p.cuts = [];
pid = 0;
p.fixedvariable = [];
p.fixdir = '';
lastUpper = upper;
oldp = p;

% Detect if pure integer or mixed
if length(p.integer_variables)+length(p.binary_variables) == length(p.c)
    p.all_integers = 1;
else
    p.all_integers = 0;
end
p.noninteger_variables = setdiff(1:length(p.c),[p.integer_variables p.binary_variables p.semicont_variables]);
poriginal.noninteger_variables = p.noninteger_variables;

% Trivial stuff in SDP cone, such as constant non-zero in row forcing
% diagonal to be non-zero
p = addImpliedSDP(p);

% Resuse some code from cutsdp to add simple cuts required for SDP
% feasibility for problems with some trivial symmetries
p = detectAndAdd3x3SymmetryGroups(p,[],0);

% Detect some more simple cuts
p = detectAndAdd3x3SDPGUBGroups(p);

% Extract more info from some combinatorial AND structure
p = cross_binary_product_cardinality(p);

% Might have detect stuff, copy back to global model
poriginal.lb = max(poriginal.lb,p.lb);
poriginal.ub = min(poriginal.ub,p.ub);

% Prepare to run...
History.feasibility = [];
History.lower = [];
History.upper = [];
History.stacksize = [];
allFeasibleSolutions = [];
allRelaxedSolutions = [];
sosgroups = [];
sosvariables = [];
unknownErrorCount = 0;
numericalProblemsCount= 0;

% Generalized upper solver format
upperSolversList = strsplit(uppersolver,',');

% A small machine for heuristics on selecting
% up- or down-node first from the stack. Without any 
% found feasible solutions, we are agnostic and alternate.
% As we gather info about where we find solutions, we 
% start picking more from one of them first.
NodeSelector.Current = 'up';
NodeSelector.Count = 1;
NodeSelector.directionWhenFoundFeasible = 0;
NodeSelector.directionWhenFoundBetter = 0;

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
    disp(['* Upper solver   : ' uppersolver]);
    disp(['* Max time       : ' num2str(p.options.bnb.maxtime)]);
    disp(['* Max iterations : ' num2str(p.options.bnb.maxiter)]);
    nB = length(p.binary_variables);
    nI = length(p.integer_variables);
    nC = length(p.c) - nB - nI;
    nKS = length(p.knapsack.a);
    nGUB = length(find(p.knapsack.type == 3));
    nCARD = length(find(p.knapsack.type == 2));
    nGEN = length(p.knapsack.type)-nGUB - nCARD;
    % Don't count the cover cardinalities derived from knapsacks
    % They are not added to model
    nGEN = nGEN - length(find(p.knapsack.type == 4));
    disp(['* ' num2str(nB) ' binaries, ' num2str(nI) ' integers, ' num2str(nC) ' continuous.']);
    disp(['* ' num2str(nKS) ' knapsacks (' num2str(nGUB) ' GUBs, ' num2str(nCARD) ' cardinality, ' num2str(nGEN) ' general)']);
    
    if possiblynonconvex & p.options.warning
        disp(' ');
        disp('Warning : The continuous relaxation may be nonconvex. This means ');
        disp('that the branching process is not guaranteed to find a');
        disp('globally optimal solution, since the lower bound can be');
        disp('invalid. Hence, do not trust the bound or the gap...')        
    end
end

if p.options.bnb.maxiter >= 0 && p.options.bnb.verbose;            disp(' Node       Upper       Gap(%)     Lower     Open  Cuts   Elapsed time');end;

% Some tricks and strategies are performed for certain objectives
p = detect_special_objectives(p);
poriginal.UnitContinuousCost = p.UnitContinuousCost;

% FIXME do when detecting gubs
p.gubs = zeros(length(p.c),1);
for i = 1:length(p.c)    
    for k = find(p.knapsack.type == 3)
        if ismember(i,p.knapsack.variables{k})
            p.gubs(i) = k;
            break
        end
    end
end

% Detect global cardinality on all binary
poriginal = extract_global_cardinality(poriginal);
p.globalcardinality = poriginal.globalcardinality;

%covercandidates = [];
%violationscount = 0;
globalcuts = [];
%cuts.knapsack = [];
%cuts.sdpknapsack = [];
%upper = 250;
sdpCuts = emptyNumericalModel;
 
while unknownErrorCount < 10 && ~isempty(node) && (etime(clock,bnbsolvertime) < p.options.bnb.maxtime) && (solved_nodes < p.options.bnb.maxiter) && (isinf(lower) || gap>p.options.bnb.gaptol)
              
    % ********************************************
    % FIXED ALONG THE TREE...
    % ********************************************
    binary_variables = p.binary_variables;
	semicont_variables  = p.semicont_variables;
      
    % ********************************************
    % ASSUME THAT WE WON'T FATHOME
    % ********************************************
    keep_digging = 1;
    Message = 'Node solved succesfully.';
    
    % We will use the change in fixed variables as
    % an indicator when selecting node,so remember 
    % where we started
    p.nfixed = nnz(p.lb == 1) + nnz(p.ub == 0);
        
    % *************************************
    % SOLVE NODE PROBLEM
    % *************************************
    if any(p.ub < p.lb - 1e-12)
        % Trivially infeasible
        x = zeros(length(p.c),1);
        output.Primal = x;
        output.problem = 1;
    else
        p.x_min = x_min;
        relaxed_p = integer_relax(p);       
        
        % Solve node relaxation         
        if p.UnitContinuousCost && sdpCuts.K.l > 0
            upper_bound = upper_from_sdpcuts(p,sdpCuts);
            if upper_bound > upper                
                output.problem = 1;
                output.Primal = [];
            else
                output = bnb_solvelower(lowersolver,addInequality(relaxed_p,globalcuts),upper,lower,x_min,allFeasibleSolutions);
            end
        else
            output = bnb_solvelower(lowersolver,addInequality(relaxed_p,globalcuts),upper,lower,x_min,allFeasibleSolutions);
        end
        
        if (output.problem == 12 || output.problem == 2) && ~(isinf(p.lower) || isnan(p.lower))
            output.problem = 1;
            Message = 'Infeasible node.'
        elseif output.problem == -1
            % This is the dreaded unknown state from mosek. Try without
            % objective to see if it is infeasible?            
            outputtest = bnb_solvelower(lowersolver,remove_objective(ptest),upper,lower,x_min,allFeasibleSolutions);
            if outputtest.problem == 1
                output.problem = 1;
                Message = 'Infeasible node.'
            else
                output.problem = -1;
            end
        end                       
        
        if output.problem == 9
            unknownErrorCount = unknownErrorCount + 1;
            Message = 'Unknown error, will try to recover.'
        end
       
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
            x  = setnonlinearvariables(p,output.Primal);
            if (p.K.l>0) && any(p.F_struc(p.K.f+1:p.K.f+p.K.l,:)*[1;x]<-1e-5)
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
    if output.problem == 0
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
    else
        % If we have numerical problems, we cannot trust current solution
        non_integer_binary = binary_variables;
        non_integer_integer = integer_variables;
        non_semivar_semivar = semicont_variables;
    end
    
    x  = setnonlinearvariables(p,x);
    
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
            
    if output.problem==0 || output.problem==3 || output.problem==4 || output.problem==5
        cost = computecost(f,c,Q,x,p);
        
        if output.problem~=1
            if output.problem == 3 || output.problem == 4 || output.problem == 5
                cost = -inf;
            end
            
            if isnan(lower)
                lower = cost;
            end
            
            if cost <= upper & ~(isempty(non_integer_binary) & isempty(non_integer_integer) & isempty(non_semivar_semivar))                              
                for k = 1:length(upperSolversList)                    
                    [upper1,x_min1,globalcuts] = feval(upperSolversList{k},p,upper,x,poriginal,output,lower,globalcuts);
                    if ~isinf(upper1)                       
                        if strcmp(p.fixdir,'up')
                            NodeSelector.directionWhenFoundFeasible = NodeSelector.directionWhenFoundFeasible + 1;
                        else
                            NodeSelector.directionWhenFoundFeasible = NodeSelector.directionWhenFoundFeasible - 1;
                        end                           
                    end
                    if upper1 < upper 
                        % Where are we finding solutions?
                        if strcmp(p.fixdir,'up')
                            NodeSelector.directionWhenFoundBetter =  NodeSelector.directionWhenFoundBetter + 1;                            
                        else
                            NodeSelector.directionWhenFoundBetter =  NodeSelector.directionWhenFoundBetter - 1;                            
                        end
                        x_min = x_min1;
                        allFeasibleSolutions = [allFeasibleSolutions x_min1];
                        upper = upper1;
                        if length(stack.nodes)>0
                            [stack,stacklower] = prune(stack,upper,p.options,solved_nodes,p,allFeasibleSolutions);
                            lower = min(lower,stacklower);
                        end
                        [p,stack] = simpleConvexDiagonalQuadraticPropagation(p,upper,stack);
                        stack = probeStackForNewUpper(stack,upper,p);
                        p = probeForNewUpper(p,upper);
                    elseif ~isinf(upper1) && upper1 == upper && norm(x_min-x_min1) > 1e-4
                        % Yet another solution with same value
                        allFeasibleSolutions = [allFeasibleSolutions x_min1];
                    end                
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
        case {-1,3,4,5,11}
            Message = 'Numerical problems, will try to recover.';
            % Solver behaved weird. Make sure we continue digging if possible          
            keep_digging = length(integer_variables)>0 || length(binary_variables)>0 || length(semicont_variables)>0;
            feasible = 1;
            cost = lower;
            x = p.lb + (p.ub-p.lb)*(1/pi);
            unbounded = find(isinf(p.ub) & isinf(p.lb));
            if ~isempty(unbounded)
                x(unbounded) = 0;
            end            
        case 0
            if can_use_ceil_lower
                lower = ceil(lower-1e-6);
            end
        case {1,12,-4,22,24}
            Message = 'Infeasible node.';
            keep_digging = 0;
            cost = inf;
            feasible = 0;
            
        case 2
            cost = -inf;
        otherwise
            % This part has to be much more robust because this could be a
            % numerical problem leading to a non-trustworthy solution
            cost = f+c'*x+x'*Q*x;
    end
    
    % **************************************
    % YAHOO! INTEGER SOLUTION FOUND
    % **************************************
    if isempty(non_integer_binary) & isempty(non_integer_integer)  & isempty(non_semivar_semivar) & ~(output.problem == -1) &  ~(output.problem == 4) & ~(output.problem == 2)        
        Message = 'Integer solution in relaxation.';
        if strcmp(p.fixdir,'up')
            NodeSelector.directionWhenFoundFeasible =  NodeSelector.directionWhenFoundFeasible + 1;
        else
            NodeSelector.directionWhenFoundFeasible =  NodeSelector.directionWhenFoundFeasible - 1;
        end
        if (cost<upper) & feasible
            if strcmp(p.fixdir,'up')
                NodeSelector.directionWhenFoundBetter =  NodeSelector.directionWhenFoundBetter + 1;
            else
                NodeSelector.directionWhenFoundBetter =  NodeSelector.directionWhenFoundBetter - 1;                
            end
            x_min = x;
            upper = cost;
            allFeasibleSolutions = [allFeasibleSolutions x_min];
            [stack,lower] = prune(stack,upper,p.options,solved_nodes,p,allFeasibleSolutions);            
            [p,stack] = simpleConvexDiagonalQuadraticPropagation(p,upper,stack);            
            % Remove trivially bad variables
            stack = probeStackForNewUpper(stack,upper,p);            
            p = probeForNewUpper(p,upper);
        end
        p = adaptivestrategy(p,upper,solved_nodes);
        keep_digging = 0;        
    elseif output.problem == 0 && ~isempty(p.binary_variables)
        
        % Not integer feasible, 
        % Try to add some GLOBAL cuts 
        % Start with knapsack
        if p.options.bnb.cut.knapsack.cover
            covers_knap1 = knapsack_add_cover_cut(p,output.Primal,'crowder',upper);
            covers_knap2 = knapsack_add_cover_cut(p,output.Primal,'gu',upper);
        else
            covers_knap1 = [];
            covers_knap2 = [];
        end
        
        % Derive knapsack from SDP cone
        if p.options.bnb.cut.sdpknapsack.cover
            allRelaxedSolutions = [allRelaxedSolutions output.Primal];
            x_trial = output.Primal;            
            [p_sdpknapsackcovers,p_rawcuts,p_sdpcuts] = create_sdpknapsack_cuts(p,x_trial,poriginal,upper);
            if p_sdpknapsackcovers.K.l > 0
                globalcuts = [globalcuts;p_sdpknapsackcovers.F_struc];
                sdpCuts = addInequality(sdpCuts,p_sdpcuts.F_struc);
            end
            if p_rawcuts.K.l>0
                % Propagate this SDP cut on the global model
                % using Glovers-Sherali second-order 
                [L,U] = propagate_second_order_cover(p_rawcuts,poriginal);
                if any(~isinf(L)) || any(~isinf(U))
                    % Apply to node and global model
                    p.lb = max(p.lb,L);
                    p.ub = min(p.ub,U); 
                    stack = tightenStackBounds(stack,L,U);                    
                    poriginal.lb = max(poriginal.lb,L);
                    poriginal.ub = min(poriginal.ub,U);              
                end
                % Do the same thing on the local model
                ptemp = p;ptemp.binary_variables = poriginal.binary_variables;
                [L,U] = propagate_second_order_cover(p_rawcuts,ptemp);
                if any(~isinf(L)) || any(~isinf(U))
                    if all(p.lb-1e-3<=x_min) & all(p.ub+1e-3>=x_min)
                        if ~(all(L-1e-3<=x_min) & all(U+1e-3>=x_min))
                            1
                        end
                    end
                    p.lb = max(p.lb,L);
                    p.ub = min(p.ub,U);                      
                end
            end
        end            
        globalcuts = unique([globalcuts;covers_knap1;covers_knap2],'rows');      
%        spy(globalcuts);drawnow
    end
    
    % **************************************
    % Stop digging if it won't give sufficient improvement anyway
    % **************************************
    if cost>upper*(1-p.options.bnb.gaptol)
        keep_digging = 0;
        if output.problem ~=1
            Message = '-> Node terminated from bound.';
        end
    end
    
    History.feasibility(end+1) = feasible;
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
        
        p0.feasible = 1;
        p1.feasible = 1;
        p0 = propagate_cardinality(p0);
        p1 = propagate_cardinality(p1);
        p0 = propagate_atleast(p0);
        p1 = propagate_atleast(p1);
        p0 = propagate_downforce(p0);
        p1 = propagate_downforce(p1);
        p0 = propagate_upforce(p0);
        p1 = propagate_upforce(p1);
        p0 = propagate_binary_product(p0);
        p1 = propagate_binary_product(p1);
        p0 = propagate_cardinality(p0);
        p1 = propagate_cardinality(p1);
            
        p0.deltafixed = nnz(p0.lb==p0.ub) - p0.nfixed;
        p1.deltafixed = nnz(p1.lb==p1.ub) - p1.nfixed;
        
        node1 = newNode(p1,globalindex,p1.fixdir,TotalIntegerInfeas,TotalBinaryInfeas,1-(x(globalindex)-floor(x(globalindex))),pid);
        pid = pid + 1;
        node0 = newNode(p0,globalindex,p0.fixdir,TotalIntegerInfeas,TotalBinaryInfeas,1-(x(globalindex)-floor(x(globalindex))),pid);
        pid = pid + 1;
                                  
        % Make sure we don't push trivially poor stuff to stack, so reuse
        % pruning code by creating temporary stacks first
        tempstack = stackCreate;
        if p1.feasible
            tempstack = push(tempstack,node1);
        end        
        if p0.feasible            
            tempstack = push(tempstack,node0);
        end
                
        tempstack = prune(tempstack,upper,p.options,solved_nodes,p,allFeasibleSolutions);
        stack = mergeStack(stack,tempstack);         
    end
        
    if stackLength(stack)>0
        lower = stackLower(stack);
        if isinf(lower) && lower > 0 && ~isinf(upper)
            lower = upper;
        else
            if can_use_ceil_lower
                lower = ceil(lower-1e-6);
            end
        end
    end
    
    % Close current and proceed to next if gap large
    gap = abs((upper-lower)/(1e-3+abs(upper)+abs(lower)));
    [node,stack,NodeSelector] = pull(stack,p.options.bnb.method,x_min,upper,NodeSelector);
    if ~isempty(node)
        p = copyNode(p,node);
    else
        % There are no nodes left
        if ~isinf(upper)
            gap = 0;
        end
    end    
    if isnan(gap)
        gap = inf;
    end    
    History.lower(end+1) = lower;
    History.upper(end+1) = upper;
    History.stacksize(end+1) = stackLength(stack);
      
    if p.options.bnb.verbose
        if mod(solved_nodes-1,p.options.print_interval)==0 || isempty(node) || (gap == 0) || (lastUpper-1e-6 > upper)
            if p.options.bnb.plot
                hold off
                subplot(1,3,1);               
                l = stairs([lowerhist' upperhist']);set(l,'linewidth',2);
                grid on
                title('Upper/lower bounds')
                subplot(1,3,2);                
                l = stairs(stacksizehist);set(l,'linewidth',2);
                grid on
                title('Open nodes')
                drawnow
                subplot(1,3,3);   
                hist(getStackLowers(stack),25);
                title('Histogram lower bounds')
                drawnow
            end
            if lastUpper > upper
                fprintf(' %4.0f : %12.3E  %7.2f   %12.3E  %2.0f   %2.0f  %8.1f    %s \n',solved_nodes,upper,100*gap,lower,stackLength(stack),size(globalcuts,1),etime(clock,bnbsolvertime),'-> Found improved solution!');
            else
                fprintf(' %4.0f : %12.3E  %7.2f   %12.3E  %2.0f   %2.0f  %8.1f    %s \n',solved_nodes,upper,100*gap,lower,stackLength(stack),size(globalcuts,1),etime(clock,bnbsolvertime),Message);
            end
        end
    end
    lastUpper = upper;    
end
if p.options.bnb.maxiter >= 0 && p.options.bnb.verbose;showprogress([num2str2(solved_nodes,3)  ' Finishing.  Cost: ' num2str(upper) ],p.options.bnb.verbose);end
if p.options.bnb.maxiter >= 0 && p.options.bnb.verbose
    if ~isempty(p.knapsack.a)
        disp(['* Added ' num2str(size(globalcuts,1)) ' global cover cuts']);
    end
end
if unknownErrorCount == 10
     diagnostics = 9;
end

% **********************************
%% BRANCH VARIABLE
% **********************************
function [index,whatsplit,globalindex] = branchvariable(x,integer_variables,binary_variables,options,x_min,Weight,p)
all_variables = [integer_variables(:);binary_variables(:)];

if ~isempty(p.sosvariables) && isempty(setdiff(all_variables,p.sosvariables)) & strcmp(options.bnb.branchrule,'sos')
    % All variables are in SOS1 constraints
    for i = 1:length(p.sosgroups)
        dist(i) = (sum(x(p.sosgroups{i}))-max(x(p.sosgroups{i})))/length(p.sosgroups{i});
    end
    % Which SOS to branch on
    [val,index] = max(dist);
    whatsplit = 'sos1';
    globalindex = index;  
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
        indic = (1+min(10,abs(p.c(all_variables)))).*(abs(x(all_variables)-round(x(all_variables))));
        [val,index] = max(indic);
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
p0.binary_variables = new_binary;

p1.ub(variable)=1;
p1.lb(variable)=1;
if length(friends) > 1
    p1.ub(friends)=0;
    p1.lb(friends)=0;
end

p0.fixdir = 'down';
p1.fixdir = 'up';
p1.binary_variables = new_binary;
p1.lower = lower;
p1.depth = p.depth+1;

% if rand<.5
%     t=p0;p0=p1;p1=t;
% end

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

p0.fixdir = 'down';
p1.fixdir = 'up';

function [p0,p1] = sos1split(p,x,index,lower,options,x_min)

v = p.sosgroups{index};
n = ceil(length(v)/2);
v1 = v(randperm(length(v),n));
v2 = setdiff(v,v1);

% In first node, set v2 to 0 and v1 to sosgroup
p0 = p;p0.lower = lower;
p0.sosgroups{index} = v1;
p0.ub(v2) = 0;

% In second node, set v1 to 0 and v1 to sosgroup
p1 = p;p1.lower = lower;
p1.sosgroups{index} = v2;
p1.ub(v1) = 0;

error('FIX down/up flag')

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

p0.fixdir = 'down';
p1.fixdir = 'up';

function s = num2str2(x,d,c)
if nargin==3
    s = num2str(x,c);
else
    s = num2str(x);
end
s = [repmat(' ',1,d-length(s)) s];


function [stack,lower] = prune(stack,upper,options,solved_nodes,p,allSolutions)
% *********************************
% PRUNE STACK W.R.T NEW UPPER BOUND
% *********************************
if stackLength(stack)>0 && ~isinf(upper)
    if length(p.binary_variables) + length(p.integer_variables) == length(p.c) && all(p.c == fix(p.c)) && nnz(p.Q)==0 && isempty(p.evalMap) && nnz(p.variabletype)==0       
        L = stack.lower;
        tooLarge = find(~isinf(L) & L>=upper-0.999);
    else     
        L = stack.lower;
        tooLarge = find(~isinf(L) & L>=upper*(1-options.bnb.prunetol));
    end
    if ~isempty(tooLarge)        
        stack.nodeCount = stack.nodeCount - length(tooLarge);
        stack.lower(tooLarge) = inf;        
    end
end

% Prune simple linear model w.r.t bound constraints only
if nnz(p.Q) == 0 && isempty(p.evalMap) && nnz(p.variabletype)==0
    tooLarge = [];
    for i = find(~isinf(stack.lower))
        pi = stack.nodes{i};
        neg = find(p.c < 0);
        pos = find(p.c > 0);
        obj = p.f + p.c(pos)'*pi.lb(pos) + p.c(neg)'*pi.ub(neg);
        if obj >= upper
            tooLarge = [tooLarge i];
        end
    end
    if ~isempty(tooLarge)
        stack.nodeCount = stack.nodeCount - length(tooLarge);
        stack.lower(tooLarge) = inf;  
    end
    if p.LinearBinaryPositiveCost && ~isinf(upper)
        % obj = sum (c_i >=0)*xi xi binary
        % if any ci is larger than upper, ub has to be 0
        % if lb is 1, it is infeasible
        r = find(p.c > upper);
        if ~isempty(r)
            for i = find(~isinf(stack.lower))
                if any(stack.nodes{i}.lb(r))
                    stack.lower(i) = inf;
                else
                    stack.nodes{i}.ub(r) = 0;
                end
            end
        end
    end
end

if stack.nodeCount > 0
    lower = min(stack.lower);
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

function p = copyNode(p,node)
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
p.pid = node.pid;
p.sosgroups = node.sosgroups;
p.sosvariables = node.sosvariables;
%p.atmost = node.atmost;
%p.localatmost = node.localatmost;

function stack = stackCreate
stack.nodes = {};
stack.lower = [];
stack.nodeCount = 0;

function stack = push(stack,p)
stack.nodes{end + 1} = p;
stack.lower(end + 1) = p.lower;
stack.nodeCount = stack.nodeCount + 1;   

function stack1 = mergeStack(stack1,stack2)
for i = 1:1:length(stack2.nodes)
    if ~(isinf(stack2.lower(i)) && stack2.lower(i)>0)
        stack1.nodes{end + 1} = stack2.nodes{i};
        stack1.lower(end + 1) = stack2.lower(i);
        stack1.nodeCount = stack1.nodeCount + 1;   
    end
end

function stack = compressStack(stack)

used = find(~(isinf(stack.lower) & stack.lower > 0));
stack.lower = stack.lower(used);
stack.nodes = {stack.nodes{used}};

function [p,stack,NodeSelector] = pull(stack,method,x_min,upper,NodeSelector)

if stackLength(stack) > 0
    if length(stack.lower)>0 && all(isinf(stack.lower))
        p = [];
        return
    end
    if numel(stack.lower) > 100 && nnz(isinf(stack.lower)) > 0.25*numel(stack.lower)
        stack = compressStack(stack);
    end
    switch method
        case {'depth','depthfirst','depthbreadth','depthproject','depthbest'}
            depths = getStackDepths(stack);
            [i,j] = max(depths);
            % Silly for backward testing compatibility. Stack order has
            % changed, to be able to compare some examples, make sure we
            % traverse the tree in exactly the same was as before on ties
            candidates = find(depths == i);
            fixed = [];
            p = [];
            for j = candidates
                %fixed = [fixed nnz(stack.nodes{j}.lb == stack.nodes{j}.ub)];
                fixed = [fixed stack.nodes{j}.deltafixed];
            end
            if ~all(fixed == fixed(1))
                [~,loc] = max(fixed);
                selected_node = candidates(loc);
                p = getStackNode(stack,selected_node);               
            else
                % Seclect based on up/down
                for j = candidates
                    if isequal(stack.nodes{j}.fixdir, NodeSelector.Current)
                        selected_node = j;
                        p = getStackNode(stack,selected_node);                       
                        NodeSelector.Count = NodeSelector.Count - 1;
                        if NodeSelector.Count == 0
                            NodeSelector.Count = 1;
                            if strcmp(NodeSelector.Current,'up')
                                NodeSelector.Current = 'down';
                            else
                                NodeSelector.Current = 'up';
                            end
                        end
                        break
                    end
                end
            end
            if isempty(p)
                selected_node = candidates(end);
                p = getStackNode(stack,selected_node);
            end                        
            %j = max(find(depths == i));
            %p = getStackNode(stack,j);
            %stack = removeStackNode(stack,j);
            stack = removeStackNode(stack,selected_node);                     
            
        case 'project'
            error
            [i,j]=min([stack.projection]);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);
            
        case 'breadth'
            error
            [i,j]=min([stack.depth]);
            p=stack(j);
            stack = stack([1:1:j-1 j+1:1:end]);
            
        case 'best'                       
            lowers = getStackLowers(stack);
            [i,j] = min(lowers);       
            candidates = find(abs(lowers-i)<=1e-10);
          
            fixed = [];
            p = [];
            for j = candidates
                fixed = [fixed nnz(stack.nodes{j}.lb == stack.nodes{j}.ub)];           
            end
            if ~all(fixed == fixed(1))
                [~,loc] = max(fixed);
                selected_node = candidates(loc);
                p = getStackNode(stack,selected_node);
            else               
                for j = candidates
                    if isequal(stack.nodes{j}.fixdir, NodeSelector.Current)
                        selected_node = j;
                        p = getStackNode(stack,selected_node);
                        NodeSelector.Count = NodeSelector.Count - 1;
                        if NodeSelector.Count == 0
                            NodeSelector.Count = 1;
                            if strcmp(NodeSelector.Current,'up')
                                NodeSelector.Current = 'down';
                                NodeSelector.Count = max(1,1 - NodeSelector.directionWhenFoundFeasible - NodeSelector.directionWhenFoundBetter);
                            else
                                NodeSelector.Current = 'up';
                                NodeSelector.Count = max(1,1 + NodeSelector.directionWhenFoundFeasible + NodeSelector.directionWhenFoundBetter);
                            end
                        end
                        break
                    end
                end
            end
            if isempty(p)                
                selected_node = candidates(end);
                p = getStackNode(stack,selected_node);               
            end    
            
            % Silly for backward testing compatibility. Stack order has
            % changed, to be able to compare some examples, make sure we
            % traverse the tree in exactly the same was as before on ties
           % j = max(find(lowers == i));
           % p = getStackNode(stack,j);
            stack = removeStackNode(stack,selected_node);
            
        otherwise
    end
else
    p = [];
end

function n = stackLength(stack)
n = stack.nodeCount;

function [L,pos,N] = stackLower(stack)
if stack.nodeCount > 0
    if nargout == 1
        L = min(stack.lower);
    elseif nargout == 2
        [L,pos] = min(stack.lower);
    elseif nargout == 3
        [L,pos] = min(stack.lower);
        N = stack.nodes(pos);
    end
else
    L = nan;
end

function D = getStackDepths(stack)
D = -inf(1,length(stack.nodes));
for i = find(~(isinf(stack.lower) & stack.lower>0))
    D(i) = stack.nodes{i}.depth;
end

function L = getStackLowers(stack)
L = stack.lower;

function N = getStackNode(stack,j)
N = stack.nodes{j};

function stack = removeStackNode(stack,j)
stack.nodeCount = stack.nodeCount - length(j);
stack.lower(j) = inf;   

function node1 = newNode(p1,globalindex,fixdir,TotalIntegerInfeas,TotalBinaryInfeas,IntInfeas,pid)
node1.lb = p1.lb;
node1.ub = p1.ub;
node1.depth = p1.depth;
node1.lower = p1.lower;
node1.fixedvariable = globalindex;
node1.fixdir = fixdir;
node1.TotalIntegerInfeas = TotalIntegerInfeas;
node1.TotalBinaryInfeas = TotalBinaryInfeas;
node1.IntInfeas = IntInfeas;
node1.x0 = p1.x0;
node1.binary_variables = p1.binary_variables;
node1.semicont_variables = p1.semicont_variables;
node1.semibounds = p1.semibounds;
node1.pid = pid;
node1.sosgroups = p1.sosgroups;
node1.sosvariables = p1.sosvariables;
%node1.atmost = p1.atmost;
%node1.localatmost = p1.localatmost;
node1.deltafixed = p1.deltafixed;
function [p,stack] = simpleConvexDiagonalQuadraticPropagation(p,upper,stack)
if ~isinf(upper) && isempty(p.evalMap) && ~any(p.c) && nnz(p.Q-diag(diag(p.Q)))==0
    d = diag(p.Q);
    if all(d)>=0
        i = find(d);
        U = sqrt(upper./d(i));
        p.lb(i) = max(p.lb(i),-U);
        p.ub(i) = min(p.ub(i),U);
        p = update_integer_bounds(p);
        fake.integer_variables = p.integer_variables;
        fake.binary_variables = p.binary_variables;
        fake.implied_integers = p.implied_integers;
        for j = 1:length(stack.nodes)                        
            stack.nodes{j}.lb(i) = max(stack.nodes{j}.lb(i),-U);
            stack.nodes{j}.ub(i) = min(stack.nodes{j}.ub(i),U);
            fake.lb =  stack.nodes{j}.lb;
            fake.ub =  stack.nodes{j}.ub;
            fake = update_integer_bounds(fake);
            stack.nodes{j}.lb = fake.lb;
            stack.nodes{j}.ub = fake.ub;
        end
    end
end

function stack = probeStackForNewUpper(stack,upper,p)
for i = 1:length(stack.nodes)
    if ~isinf(stack.lower(i))
        stack.nodes{i} = probeForNewUpper(stack.nodes{i},upper,p);
    end
end

function s = probeForNewUpper(s,upper,p)
if nargin == 2
    p = s;
end
% FIXME generalize for negative etc
if p.LinearBinaryPositiveCost
    % Compute possible indivual contribution to cost
    % Only tests those that can be 0 or 1
    probe = find(s.lb(s.binary_variables) < s.ub(s.binary_variables));
    probe = s.binary_variables(probe);
    % This is as low we can go if we keep everything at lb
    trivialBound = sum(p.c.*s.lb);
    % Check what happens if we fix a variable to 1
    must_be_zero = [];
    for j = probe
        % Fixing leads to poor bound
        if trivialBound + p.c(j) > upper
            must_be_zero = [must_be_zero j];
        end
    end
    if ~isempty(must_be_zero)
        s.ub(must_be_zero) = 0;
    end
end

function stack = tightenStackBounds(stack,L,U)
for i = 1:length(stack.nodes)
    if ~isinf(stack.lower(i))
        stack.nodes{i}.lb = max(stack.nodes{i}.lb,L);
        stack.nodes{i}.ub = min(stack.nodes{i}.ub,U);
        if any(stack.nodes{i}.lb > stack.nodes{i}.ub)
           stack.lower(i) = inf;
        end
    end
end