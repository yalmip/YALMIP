function [x_min,solved_nodes,lower,upper,profile,diagnostics] = bnb_branch_and_bound(p)

% Prepare root node
p.options.savesolveroutput = 0;
p.options.savesolverinput = 1; % We save this to re-use some stuff in fmincon
p.options.saveduals = 0;
p.options.dimacs = 0;
p.getsolvertime = p.options.bnb.profile;
p.options.allowsmashing = 1;
p.depth = 0;
p.lower = NaN;
p.cuts = [];
p.fixdir = '';

diagnostics = 0;
bnbsolvertime = clock;
profile.local_solver_time = 0;
lower = NaN;
solved_nodes = 0;
gap = inf;
node = 1;
stack = stackCreate;
lowersolver = p.solver.lower.call;
uppersolver = p.options.bnb.uppersolver;
History.feasibility = [];
History.lower = [];
History.upper = [];
History.stacksize = [];
allFeasibleSolutions = [];
allRelaxedSolutions = [];
unknownErrorCount = 0;
numericalProblemsCount= 0;
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

% Detect common case with 1 continuous variable in SDP only
p = addSDPextendable(p);

% p = detectMonotoneSDPObjective(p);
p = detect_knapsack(p);
p = presolve_cliquestrengthen(p);

[p,x_min,upper] = initializesolution(p,[],inf);
if ~isinf(upper) && p.options.bnb.verbose
    disp(['* Heuristics found solution with objective ' num2str(upper)]);
end

pid = 0;
lastUpper = upper;
oldp = p;

% Remember if pure integer or mixed
p.all_integers  = (length(p.integral_variables) == length(p.c));

% Resuse some code from cutsdp to add simple cuts required for SDP
% feasibility for problems with some trivial symmetries
p = detectAndAdd3x3SymmetryGroups(p,[],0);

% Detect some more simple cuts
p = detectAndAdd3x3SDPGUBGroups(p);

% Extract more info from some combinatorial AND structure
p = cross_binary_product_cardinality(p);

p = lpcardinality(p);

if p.options.bnb.verbose
    
    pc = p.problemclass;
    non_convex_obj = pc.objective.quadratic.nonconvex | pc.objective.polynomial;
    non_convex_constraint =  pc.constraint.equalities.quadratic | pc.constraint.inequalities.elementwise.quadratic.nonconvex;
    non_convex_constraint =  non_convex_constraint | pc.constraint.equalities.polynomial | pc.constraint.inequalities.elementwise.polynomial;
    
    possiblynonconvex = non_convex_obj | non_convex_constraint;
    if ~isequal(p.solver.lower.version,'')
        p.solver.lower.tag = [p.solver.lower.tag '-' p.solver.lower.version];
    end
    
    % Number of effective binaries is not same as those
    % kept in list p.binary_variables, as some are presolved
    % (and marked as such in removableVariables)
    % away and thus not really in the problem any longer
    nB = length(p.binary_variables);
    nI = length(p.integer_variables);
    nC = length(p.noninteger_variables);
    nKS = length(p.knapsack.a);
    nGUB = length(find(p.knapsack.type == 3));
    nCARD = length(find(p.knapsack.type == 2));
    nGEN = length(p.knapsack.type)-nGUB - nCARD;
    % Don't count the cover cardinalities derived from knapsacks
    % They are not added to model
    nGEN = nGEN - length(find(p.knapsack.type == 4));
    nIMPL = length(intersect(p.binary_variables,p.implied_integers));
    if nIMPL > 0
        disp(['* Presolved model: ' num2str(nB) ' binaries (' num2str(nIMPL) ' implied), ' num2str(nI) ' integers, ' num2str(nC) ' continuous.']);
    else
        disp(['* Presolved model: ' num2str(nB) ' binaries, ' num2str(nI) ' integers, ' num2str(nC) ' continuous.']);
    end
    disp(['*                  ' num2str(p.K.f+p.K.l) ' linear rows, ' num2str(sum(p.K.q)+sum(p.K.s.^2)) ' conic rows']);
    
    disp(['* ' num2str(nKS) ' knapsacks (' num2str(nGUB) ' GUBs, ' num2str(nCARD) ' cardinality, ' num2str(nGEN) ' general)']);
    
    if possiblynonconvex && p.options.warning
        disp(' ');
        disp('Warning : The continuous relaxation may be nonconvex. This means ');
        disp('that the branching process is not guaranteed to find a');
        disp('globally optimal solution, since the lower bound can be');
        disp('invalid. Hence, do not trust the bound or the gap...')
    end
end

% Detect global cardinality on all binary from knapsacks
p = extract_global_cardinality(p);

sdpCuts = emptyNumericalModel;
AddedCoverCuts = 0;
% We have a global pack of cuts, and then we have keep track which
% we have added to the global pool. However, we also keep track of
% which we have deemed redundant in a node (and thus deeper in the tree
% from there)
globalcuts = [];
p.cutactivate  = [];% This is global information
p.cutredundant = [];% Local node copied to children

% Some simple trial solutions
[x_min,upper] = root_node_heuristics(p,x_min,upper);
if ~isinf(upper)
    if p.options.bnb.verbose
        disp(['* Heuristics found solution with objective ' num2str(upper)]);
    end
    % Might be able to tighten some stuff based on this solution
    p = probeForNewUpper(p,upper);
end

if p.options.bnb.verbose
    disp('Node      Incumbent      Gap        Bound   Open  Cut     Time  Diagnostics');
end

% Save the root node to keep data which does not have to
% be copied to every node, and to keep global information
poriginal = p;


while unknownErrorCount < 10 && ~isempty(node) && (etime(clock,bnbsolvertime) < p.options.bnb.maxtime) && (solved_nodes < p.options.bnb.maxiter) && (isinf(lower) || gap>p.options.bnb.gaptol)
    
    % ********************************************
    % ASSUME THAT WE WON'T FATHOME
    % ********************************************
    keep_digging = 1;
    p.feasible = 1;
    Message = 'Node solved succesfully.';
    % This should be perhaps be done directly on the stack
    p = updateCardinalityFromUpper(p,upper);
    
    p_start = p;
    
    p = cut_process(p,globalcuts);
    
    % We will use the change in fixed variables as
    % an indicator when selecting node,so remember
    % where we started
    p.nfixed = nnz(p.lb == 1) + nnz(p.ub == 0);
    
    % *************************************
    % SOLVE NODE PROBLEM
    % *************************************
    if any(p.ub < p.lb - 1e-12) || p.feasible == 0
        % Trivially infeasible
        x = zeros(length(p.c),1);
        output.Primal = x;
        output.problem = 1;
        p.feasible = 0;
    else
        p.x_min = x_min;
        relaxed_p = integer_relax(p);
                
        p = checkConsistentSDPcutUpper(p,sdpCuts,upper);
        if ~p.feasible
            output.problem = 1;
            output.Primal = [];
        end
        
        % Solve node relaxation
        if p.feasible
            use = ~p.cutredundant(:) & p.cutactivate(:);
            p_node = addInequality(relaxed_p,globalcuts,use);           
            output = bnb_solvelower(lowersolver,p_node,upper,lower,x_min,allFeasibleSolutions);            
        end
        
        if (output.problem == 12 || output.problem == 2) && ~(isinf(p.lower) || isnan(p.lower))
            output.problem = 1;
            Message = 'Infeasible node.';
        elseif output.problem == -1
            % This is the dreaded unknown state from mosek. Try without
            % objective to see if it is infeasible?
            outputtest = bnb_solvelower(lowersolver,remove_objective(relaxed_p),upper,lower,x_min,allFeasibleSolutions);
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
            % FIXME experimental stuff with eliminated binaries
            % These are kept, but not in any constraint
            %x(p.removableVariables) = 0;
            x = setnonlinearvariables(p,output.Primal);
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
        non_integer_binary = abs(x(p.binary_variables)-round(x(p.binary_variables)))>p.options.bnb.inttol;
        non_integer_integer = abs(x(p.integer_variables)-round(x(p.integer_variables)))>p.options.bnb.inttol;
        if p.options.bnb.round
            x(p.binary_variables(~non_integer_binary))   = round(x(p.binary_variables(~non_integer_binary)));
            x(p.integer_variables(~non_integer_integer)) = round(x(p.integer_variables(~non_integer_integer)));
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
        non_integer_binary = p.binary_variables;
        non_integer_integer = p.integer_variables;
        non_semivar_semivar = p.semicont_variables;
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
        cost = computecost(p.f,p.c,p.Q,x,p);
        
        if output.problem~=1
            if output.problem == 3 || output.problem == 4 || output.problem == 5
                if nnz(p.Q)==0 && isempty(p.evalMap)
                    cost = sum(p.c(p.c>0).*p.lb(p.c>0))+sum(p.c(p.c<0).*p.ub(p.c<0));
                else
                    cost = -inf;
                end
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
                        poriginal = probeForNewUpper(poriginal,upper);
                        poriginal = updateCardinalityFromUpper(poriginal,upper);
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
            % keep_digging = length(integer_variables)>0 || length(binary_variables)>0 || length(semicont_variables)>0;
            keep_digging = any(p.lb(p.integral_variables)~=p.ub(p.integral_variables));
            feasible = 1;
            cost = lower;
            x = p.lb + (p.ub-p.lb)*(1/pi);
            unbounded = find(isinf(p.ub) & isinf(p.lb));
            if ~isempty(unbounded)
                x(unbounded) = 0;
            end
        case 0
            if p.IntegerCost
                lower = ceil(lower-1e-5);
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
            cost = p.f+p.c'*x+x'*p.Q*x;
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
            % Remove trivially bad nodes
            stack = probeStackForNewUpper(stack,upper,p);
            poriginal = probeForNewUpper(poriginal,upper);
            p = probeForNewUpper(p,upper);
            poriginal = updateCardinalityFromUpper(poriginal,upper);
            p = updateCardinalityFromUpper(p,upper);
        end
        p = adaptivestrategy(p,upper,solved_nodes);
        keep_digging = 0;
    elseif output.problem == 0 && ~isempty(p.binary_variables)
        
        % Not integer feasible,
        % Try to add some GLOBAL cuts
        % Start with knapsack
        new_cuts = [];
        if p.options.bnb.cuts && p.options.bnb.cut.knapsack.cover
            covers_knap1 = knapsack_add_cover_cut(p,output.Primal,'crowder',upper);
            covers_knap2 = knapsack_add_cover_cut(p,output.Primal,'gu',upper);
            new_cuts = [covers_knap1;covers_knap1];
        end
        
        % Derive some cut from SDP cone and use to propagate
        if any(p.K.s) && length(p.binary_variables)>= 0 && p.options.bnb.cuts
            % p_sdpknapsack : Result from cover analysis
            % p_binarycut   : Original SDP cut with continuous fixed
            % p_sdpcuts     : Original SDP cut
            [p_sdpknapsackcovers,p_binarycut,p_sdpcuts] = create_sdpknapsack_cuts(p,output.Primal,poriginal,upper);
            if p.options.bnb.cut.sdpknapsack.cover && p_sdpknapsackcovers.K.l > 0
                new_cuts = [new_cuts;p_sdpknapsackcovers.F_struc];
            end
            if p_binarycut.K.l>0
                sdpCuts = addInequality(sdpCuts,p_sdpcuts.F_struc);
                
                % Propgate implications from this cut on cardinality
                % The cut is valid globally
                poriginal = updateLowerCardinalityFromCut(poriginal,p_binarycut);
                poriginal = updateUpperCardinalityFromCut(poriginal,p_binarycut);
                p = propagate_global_to_node(poriginal,p);
                
                % Propagate this cut on the global model
                % using Glovers-Sherali second-order covers
                % Does not reduce #nodes but is efficient at detecting 0s
                % which leads to quicker solves in the nodes
                poriginal = iterate_propagate_second_order_cover(p_binarycut,poriginal);
                poriginal = updateLowerCardinalityFromCut(poriginal,p_binarycut);
                poriginal = updateUpperCardinalityFromCut(poriginal,p_binarycut);
                p = propagate_global_to_node(poriginal,p);
                
                % Do the same thing on the local model
                p = iterate_propagate_second_order_cover(p_binarycut,p);
                p = updateLowerCardinalityFromCut(p,p_binarycut);
                p = updateUpperCardinalityFromCut(p,p_binarycut);
                
                % Exploit binary product models (remove)
                [p,poriginal] = propagate_binary_product_on_cut(p,poriginal,p_binarycut);
            end
        end
        if ~isempty(globalcuts) && ~isempty(new_cuts)
            new_hash = new_cuts*[1;p.hash];
            old_hash = globalcuts*[1;p.hash];
            s = find(~ismember(new_hash,old_hash));
            new_cuts = new_cuts(s,:);
        end
        globalcuts = [globalcuts;new_cuts];
    end
    
    % **************************************
    % Stop digging if it won't give sufficient improvement anyway
    % **************************************
    if cost>upper*(1-p.options.bnb.gaptol) || p.feasible == 0
        keep_digging = 0;
        if output.problem ~=1
            Message = '-> Node terminated from bound.';
        end
    end
    
    History.feasibility(end+1) = feasible;
    % **********************************
    % CONTINUE SPLITTING?
    % **********************************
    if keep_digging && (cost<upper)
        
        if solved_nodes == 1
            RootNodeInfeas =  TotalIntegerInfeas+TotalBinaryInfeas;
            RootNodeCost = cost;
        end
        
        % **********************************
        % BRANCH VARIABLE
        % **********************************
        [index,whatsplit,globalindex] = branchvariable(x,p.integer_variables,p.binary_variables,p.options,x_min,[],p);
        
        p = activateCuts(p,globalcuts,x);
              
        % **********************************
        % CREATE NEW PROBLEMS
        % **********************************
        switch whatsplit
            case 'binary'
                [p0,p1,index] = binarysplit(p,x,index,cost,[]);
                
            case 'integer'
                [p0,p1] = integersplit(p,x,index,cost,x_min);
                
            case 'semi'
                [p0,p1] = semisplit(p,x,index,cost,x_min);
                
            otherwise
        end
        
        p0.feasible = 1;
        p1.feasible = 1;
        p0 = propagate_multiple(p0);
        p1 = propagate_multiple(p1);
        p0 = checkConsistentSDPcutUpper(p0,sdpCuts,upper);
        p1 = checkConsistentSDPcutUpper(p1,sdpCuts,upper);
        
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
            if p.IntegerCost
                lower = ceil(lower-1e-5);
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
    History.lower = [History.lower  lower];
    History.upper = [History.upper upper];
    History.stacksize(end+1) = stackLength(stack);
    
    if p.options.bnb.verbose
        if mod(solved_nodes-1,p.options.print_interval)==0 || isempty(node) || (gap == 0) || (lastUpper-1e-6 > upper)
            if p.options.bnb.plot
                hold off
                subplot(1,3,1);
                l = stairs([History.lower' History.upper']);set(l,'linewidth',2);
                grid on
                title('Upper/lower bounds')
                subplot(1,3,2);
                l = stairs(History.stacksize);set(l,'linewidth',2);
                grid on
                title('Open nodes')
                drawnow
                subplot(1,3,3);
                hist(getStackLowers(stack),25);
                title('Histogram lower bounds')
                drawnow
            end
            upper = full(upper);
            lower = full(lower);
            
            binary_diagnostics = createBinaryDiagnosticsString(p_start);
            
            if length(p.cutactivate)<size(globalcuts,1)
                p.cutactivate(size(globalcuts,1)) = 0;
            end
            if length(p.cutredundant)<size(globalcuts,1)
                p.cutredundant(size(globalcuts,1)) = 0;
            end
            
            if lastUpper > upper
                fprintf(' %4.0f : %12.3E  %6.2f%%  %12.3E  %2.0f  %2.0f/%2.0f %6.1f  %s  %s \n',solved_nodes,upper,100*gap,lower,stackLength(stack),nnz(~p.cutredundant & p.cutactivate),nnz(p.cutactivate),etime(clock,bnbsolvertime),binary_diagnostics,'-> Found improved solution!');
            else
                fprintf(' %4.0f : %12.3E  %6.2f%%  %12.3E  %2.0f  %2.0f/%2.0f %6.1f  %s  %s \n',solved_nodes,upper,100*gap,lower,stackLength(stack),nnz(~p.cutredundant & p.cutactivate),nnz(p.cutactivate),etime(clock,bnbsolvertime),binary_diagnostics,Message);
            end
        end
    end
    lastUpper = upper;
end
if p.options.bnb.maxiter >= 0 && p.options.bnb.verbose;showprogress([num2str2(solved_nodes,3)  ' Finishing.  Cost: ' num2str(upper) ],p.options.bnb.verbose);end
if p.options.bnb.maxiter >= 0 && p.options.bnb.verbose
    if AddedCoverCuts > 0
        disp(['* Added ' num2str(AddedCoverCuts) ' global cover cuts']);
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
        % Merit brancing on variables in GUBs
        % indic = (1+min(10,abs(p.c(all_variables)))).*(abs(x(all_variables)-round(x(all_variables))));
        if ~isempty(p.cliques.table)
            part_of_clique = any(p.cliques.table(all_variables,:),2);
            indic = (1+part_of_clique(:) + min(10,abs(p.c(all_variables)))).*(abs(x(all_variables)-round(x(all_variables))));
        else
            indic = (1+min(10,abs(p.c(all_variables)))).*(abs(x(all_variables)-round(x(all_variables))));
        end
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
function [p0,p1,variable] = binarysplit(p,x,index,lower,options)
p0 = p;
p1 = p;

variable = p.binary_variables(index);

p0.ub(variable) = 0;
p0.lb(variable) = 0;
p0.lower = lower;
p0.depth = p.depth+1;
p0.fixdir = 'down';
p0.IntInfeas = x(variable)-p0.ub(variable);

p1.ub(variable)=1;
p1.lb(variable)=1;
p1.lower = lower;
p1.depth = p.depth+1;
p1.fixdir = 'up';
p1.IntInfeas = p1.lb(variable)-x(variable);

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
p0.IntInfeas = x(variable)-p0.ub(variable);
p1.IntInfeas = p1.lb(variable)-x(variable);

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

p0.IntInfeas = 1;
p1.IntInfeas = 1;

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
    if length(p.binary_variables) + length(p.integer_variables) == length(p.c) && all(p.c == fix(p.c)) && nnz(p.Q)==0 && isempty(p.evalMap) && nnz(p.variabletype)==0 && isequal(p.K.m,0)
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
if nnz(p.Q) == 0 && isempty(p.evalMap) && nnz(p.variabletype)==0 && isequal(p.K.m,0)
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
p.cutredundant = node.cutredundant;
p.binarycardinality = node.binarycardinality;
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

function node1 = newNode(p1,globalindex,fixdir,TotalIntegerInfeas,TotalBinaryInfeas,notused,pid)
% Dynamic fields saved to stack elements
node1.lb = p1.lb;
node1.ub = p1.ub;
node1.depth = p1.depth;
node1.lower = p1.lower;
node1.binarycardinality = p1.binarycardinality;
node1.cutredundant = p1.cutredundant;
node1.fixedvariable = globalindex;
node1.fixdir = fixdir;
node1.TotalIntegerInfeas = TotalIntegerInfeas;
node1.TotalBinaryInfeas = TotalBinaryInfeas;
node1.IntInfeas = p1.IntInfeas;
node1.x0 = p1.x0;
node1.binary_variables = p1.binary_variables;
node1.semicont_variables = p1.semicont_variables;
node1.semibounds = p1.semibounds;
node1.pid = pid;
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
if nargin == 2 || isinf(upper)
    p = s;
end
if p.LinearBinaryPositiveCost && ~isinf(upper)
    % Compute possible indivual contribution to cost
    % Only tests those that can be 0 or 1
    probe = find(s.lb(s.binary_variables) < s.ub(s.binary_variables));
    probe = s.binary_variables(probe);
    % This is as low we can go if we keep everything at lb
    trivialBound = sum(p.c.*s.lb);
    % Check what happens if we fix a variable to 1
    must_be_zero = [];
    for j = probe
        % Fixing leads to poor bound (not better than current)
        if trivialBound + p.c(j) >= upper
            must_be_zero = [must_be_zero j];
        end
    end
    if ~isempty(must_be_zero)
        s.ub(must_be_zero) = 0;
    end
elseif p.LinearBinaryNegativeCost && ~isinf(upper)
    % Compute possible indivual contribution to cost
    % Only tests those that can be 0 or 1
    probe = find(s.lb(s.binary_variables) < s.ub(s.binary_variables));
    probe = s.binary_variables(probe);
    % This is as low we can go if we keep everything at ub
    trivialBound = sum(p.c.*s.ub);
    % Check what happens if we fix a variable to 0
    must_be_one = [];
    for j = probe
        % Fixing leads to poor bound (not better than current)
        if trivialBound - p.c(j) >= upper
            must_be_one = [must_be_one j];
        end
    end
    if ~isempty(must_be_one)
        s.lb(must_be_one) = 1;
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

function p = updateCardinalityFromUpper(p,upper)

if p.LinearBinaryPositiveCost && ~isinf(upper)
    fixed_binary = find(p.lb==1 & p.isbinary(:));
    free_binary = find(p.lb==0 & p.ub==1 & p.isbinary(:));
    if ~isempty(free_binary)
        c = sort(p.c(free_binary),'ascend');
        k = max(find(sum(p.c(fixed_binary))+cumsum(c) < upper));
        if isempty(k)
            if sum(p.c(fixed_binary)) > upper
                p.feasible=0;
                return
            else
                k=0;
            end
        end
    else
        k=0;
    end
    p.binarycardinality.up = min(length(fixed_binary)+k,p.binarycardinality.up);
    
elseif p.LinearBinaryCost && ~isinf(upper)
    % How many do we need to set to improve upper?
    % Find those we actually can set
    c = p.c(p.ub == 1);
    c = c(find(c));
    k = min(find(cumsum(sort(c,'ascend')) < upper));
    p.binarycardinality.down = max(k,p.binarycardinality.down);
end


function p = updateLowerCardinalityFromCut(p,p_cut)
% 5*x1+4*x2+3*x3+2*x4 >= 9 so at least 2 must be set
if p.feasible && p.binarycardinality.down~=p.binarycardinality.up
    for i = 1:p_cut.K.l
        if p_cut.F_struc(i,1) < 0
            s = p_cut.F_struc(i,1 + p.binary_variables);
            s = s(p.ub(p.binary_variables)==1);
            k = min(find(cumsum(sort(s,'descend')) >= -p_cut.F_struc(i,1)));
            if isempty(k)
                p.feasible = 0;
            else
                p.binarycardinality.down = max(p.binarycardinality.down,k);
            end
        end
    end
end
function p = updateUpperCardinalityFromCut(p,p_cut)
if p.feasible && p.binarycardinality.down~=p.binarycardinality.up
    for i = 1:p_cut.K.l
        s = p_cut.F_struc(i,1 + p.binary_variables);
        k = max(find(cumsum(sort(s,'descend')) >= -p_cut.F_struc(i,1)));
        p.binarycardinality.up = min(p.binarycardinality.up,k);
    end
end

function p = check_cliqueviolation(p)
if p.feasible
    for i = size(p.cliques.table,2):-1:1
        vars = find(p.cliques.table(:,i));
        switch nnz(p.lb(vars))
            case 0
            case 1
            otherwise
                p.feasible = 0;
                break
        end
    end
end

function p = checkConsistentSDPcutUpper(p,sdpCuts,upper)
if p.feasible && ~isinf(upper) && p.UnitContinuousCost && sdpCuts.K.l > 0
    upper_lowerbound = upper_from_sdpcuts(p,sdpCuts,upper);
    if upper_lowerbound > upper
        p.feasible = 0;        
    end
end

function p = iterate_propagate_second_order_cover(p_binarycut,p)
n_fixed = nnz(p.lb(p.binary_variables)==p.ub(p.binary_variables));
goon = 1;
while goon && p.options.bnb.cardinalitypropagations && p.feasible
    [L,U] = propagate_second_order_cover(p_binarycut,p);
    if any(L>p.ub) || any(U < p.lb)
        p.feasible = 0;
    else
        if any(~isinf(L)) || any(~isinf(U))
            p.lb = max(p.lb,L);
            p.ub = min(p.ub,U);
        end
        new_n_fixed = nnz(p.lb(p.binary_variables)==p.ub(p.binary_variables));
        if new_n_fixed == n_fixed
            goon = 0;
        end
        n_fixed = new_n_fixed;
    end
end

function binary_diagnostics = createBinaryDiagnosticsString(p)
fixed_1 = nnz(p.lb(p.binary_variables)==1);
fixed_0 = nnz(p.ub(p.binary_variables)==0);
binary_diagnostics = ['(' num2str(fixed_0) '/' num2str(fixed_1) '/' num2str(p.binarycardinality.down) '/' num2str(p.binarycardinality.up) ')'];

function p = propagate_multiple(p)
p.binarycardinality.down = max(p.binarycardinality.down,nnz(p.lb(p.binary_variables)==1));
p.binarycardinality.up = min(p.binarycardinality.up,nnz(p.ub(p.binary_variables)==1));
p = propagate_cardinality(p);
p = propagate_setcover(p);
p = propagate_downforce(p);
p = propagate_upforce(p);
p = propagate_binary_product(p);
p = propagate_cardinality(p);
p = check_cliqueviolation(p);

function p = propagate_global_to_node(poriginal,p)
p.binarycardinality.down = max(poriginal.binarycardinality.down,p.binarycardinality.down);
p.binarycardinality.up   = min(poriginal.binarycardinality.up,p.binarycardinality.up);
p.lb = max(p.lb,poriginal.lb);
p.ub = min(p.ub,poriginal.ub);
if any(p.ub < p.lb)
    p.feasible = 0;
end

function [p,poriginal] = propagate_binary_product_on_cut(p,poriginal,p_binarycut)
if p.feasible && ~isempty(p.binaryProduct)
    for i = 1:p_binarycut.K.l
        row = p_binarycut.F_struc(i,:);
        if all(row(2:end) <= 0)
            % we have c - sum xi*(ai>0) >=0
            r = row(2:end);
            e = p.binaryProduct(:,2:3);
            e = unique(e);
            if all(ismember(find(r),e))
                r = sort(r(r<0),'descend');
                n = max(find(row(1) + cumsum(r) >= 0));
                % Not all of them where in knapsack?
                n = n + length(setdiff(e,find(r)));
                % All active variables in the model is
                % contained in the z = xi*xj structures via x
                if all(ismember(find(poriginal.ub == 1),p.binaryProduct))
                    % Since at most n of the x variables can be
                    % non-zero, at most nchoosek(n,2) of the z can be
                    % non-zero
                    poriginal.binarycardinality.up = min(poriginal.binarycardinality.up,nchoosek(n,2));
                    p.binarycardinality.up = poriginal.binarycardinality.up;
                end
            end
        end
    end
end

function p = cut_process(p,globalcuts)
m = size(globalcuts,1);
if length(p.cutactivate) < m
    p.cutactivate(m) = 0;
end
if length(p.cutredundant) < m
    p.cutredundant(m) = 0;
end
if ~isempty(globalcuts)
    s = find(p.cutactivate & ~p.cutredundant);
    A = globalcuts(s,2:end);
    b = globalcuts(s,1);
    locally_redundant = find((b + [A.*(A>0)]*p.lb + [A.*(A<0)]*p.ub) >= 0);
    p.cutredundant(s(locally_redundant)) = 1;
end

function p = activateCuts(p,globalcuts,x)
if size(globalcuts,1) > 0
    % Implementing cutting plane management and selection techniques
    % Franz Wesselmann, Uwe H. Suhl
    if length(p.cutactivate)<size(globalcuts,1)
        p.cutactivate(size(globalcuts,1)) = 0;
    end
    if length(p.cutredundant)<size(globalcuts,1)
        p.cutredundant(size(globalcuts,1)) = 0;
    end
    for i = 1:size(globalcuts,1)
        a0 = globalcuts(i,1);
        a  = -globalcuts(i,2:end);
        eff(i,1) = (a*x - a0)/norm(a);
        isp(i,1) = nnz(p.c(p.binary_variables))/nnz(p.c);
        obp(i,1) = abs(a*p.c/(norm(a)*norm(p.c)));
    end
    meas = [eff isp obp]*[1/3;1/3;1/3];
    meas(eff<0) = -inf;
    meas(obp<0.1) = -inf;
    meas(find(p.cutactivate))=-inf;
    for loops = 1:4
        if any(~isinf(meas))
            [~,j] = max(meas);
            p.cutactivate(j) = 1;
            meas(j) = -inf;
            a = globalcuts(j,2:end);
            for k = 1:size(globalcuts,1)
                b = globalcuts(k,2:end);
                if (abs(a*b')/(norm(a)*norm(b)))>0.9
                    meas(k) = -inf;
                end
            end
        end
    end
end

function p = lpcardinality(p)
% Waste two LPs on figuring out binary cardinality
if any(p.isbinary)
    % Create a model with only LP stuff
    pLP=p;
    pLP.Q=pLP.Q*0;pLP.evalMap = [];
    pLP.F_struc = pLP.F_struc(1:p.K.f+p.K.l,:);
    pLP.variabletype = spalloc(1,length(pLP.c),0);
    pLP.monomtable = eye(length(pLP.c));
    pLP.K.s=0;pLP.K.e=0;pLP.K.q=[];pLP.K.p = [];
    pLP.c=pLP.c*0;
    
    try
        % We put this in a catch, in case some weird lower solver is used
        % such as a GP solver or something
        
        % Minimize sum of all binary variables
        pLP.c(p.binary_variables)=1;
        output=feval(pLP.solver.lower.call,integer_relax(pLP));
        p.binarycardinality.down = max(p.binarycardinality.down,ceil(output.Primal'*pLP.c-1e-5));
        
        % Maximize sum of all binaries
        pLP.c=pLP.c*0;pLP.c(p.binary_variables)=-1;
        output=feval(pLP.solver.lower.call,integer_relax(pLP));
        p.binarycardinality.up = min(p.binarycardinality.up, floor(-output.Primal'*pLP.c+1e-5));
    catch
    end
end