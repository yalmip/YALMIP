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

if nnz(Q)==0 && nnz(c)==1 && ~any(p.K.m)
    p.simplecost = 1;
else
    p.simplecost = 0;
end

p = detectSOS(p);
p = detectAtMost(p);
p = presolveSOS(p);
p = smashFixed(p);
p = propagate_bounds_from_qualities(p);

p.options.allowsmashing = 1;
poriginal = p;
p.cuts = [];
pid = 0;
lowerhist = [];
upperhist = [];
stacksizehist = [];
p.fixedvariable = [];
p.fixdir = '';
lastUpper = upper;
oldp = p;

if length(p.integer_variables) == length(p.c)
    p.all_integers = 1;
else
    p.all_integers = 0;
end
p.noninteger_variables = setdiff(1:length(p.c),[p.integer_variables p.binary_variables p.semicont_variables]);
poriginal.noninteger_variables = p.noninteger_variables;

p = addImpliedSDP(p);

% Resuse some code from cutsdp to add simple cuts required for SDP
% feasibility for problems with some trivial symmetries
% TODO: Clean up, refactor, generalize
if ~any(p.K.f) % Still to lazy to fix last insertion
    top = startofSDPCone(p.K);
    p.F_struc = p.F_struc';
    for i = 1:length(p.K.s)
        p.semidefinite{i}.F_struc = p.F_struc(:,top:top+p.K.s(i)^2-1)';
        top = top + p.K.s(i)^2;
    end
    p.F_struc = p.F_struc';
    p.sdpsymmetry = [];
    p = detect3x3SymmetryGroups(p);
    pp = p;pp.F_struc = [];pp.K.l = 0;pp.K.f = 0;pp.K.s = 0;
    [p,pp] = addSymmetryCuts(p,pp);
    p = addInequality(p,pp.F_struc);
    p.semidefinite=[];
end

feasibilityHistory = [];
% Save of all optimal solutions
allSolutions = [];
sosgroups = [];
sosvariables = [];
unknownErrorCount = 0;
numericalProblemsCount= 0;

% Generalized upper solver format
upperSolversList = strsplit(uppersolver,',');

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
    
    if possiblynonconvex & p.options.warning
        disp(' ');
        disp('Warning : The continuous relaxation may be nonconvex. This means ');
        disp('that the branching process is not guaranteed to find a');
        disp('globally optimal solution, since the lower bound can be');
        disp('invalid. Hence, do not trust the bound or the gap...')
        
    end
end
if p.options.bnb.verbose;            disp(' Node       Upper       Gap(%)     Lower     Open   Elapsed time');end;

while unknownErrorCount < 10 && ~isempty(node) && (etime(clock,bnbsolvertime) < p.options.bnb.maxtime) && (solved_nodes < p.options.bnb.maxiter) && (isinf(lower) || gap>p.options.bnb.gaptol)

%     for i = 1:length(p.binaryProduct)
%         p.lb(p.binaryProduct{i}.y) = max(p.lb(p.binaryProduct{i}.y),min(p.lb(p.binaryProduct{i}.x)));
%         p.ub(p.binaryProduct{i}.y) = min(p.ub(p.binaryProduct{i}.y),max(p.ub(p.binaryProduct{i}.x)));
%     end
    
    % ********************************************
    % BINARY VARIABLES ARE FIXED ALONG THE PROCESS
    % ********************************************
    binary_variables = p.binary_variables;

    % ********************************************
    % SO ARE SEMI VARIABLES
    % ********************************************
    semicont_variables  = p.semicont_variables;
    
    % ********************************************
    % ASSUME THAT WE WON'T FATHOME
    % ********************************************
    keep_digging = 1;
    message = '';
     
%     for i = 1:length(p.atmost.groups)
%         j = p.atmost.groups{i};       
%         if any(p.lb(j)==1)
%             fixed = j(find(p.lb(j)==1));
%             zerov = setdiff(j,fixed);
%             p.ub(zerov) = 0;
%         end        
%     end
     
    % *************************************
    % SOLVE NODE PROBLEM
    % *************************************
    if any(p.ub<p.lb - 1e-12)
        x = zeros(length(p.c),1);
        output.Primal = x;
        output.problem=1;
    else
        p.x_min = x_min;
        relaxed_p = p;
        relaxed_p.integer_variables = [];
        relaxed_p.binary_variables = [];
        relaxed_p.semicont_variables = [];
        relaxed_p.ub(p.ub<p.lb) = relaxed_p.lb(p.ub<p.lb);
        
        % Solve node relaxation     
        output = bnb_solvelower(lowersolver,relaxed_p,upper,lower,x_min,allSolutions);
        if (output.problem == 12 || output.problem == 2) && ~(isinf(p.lower) || isnan(p.lower))
            output.problem = 1;
        elseif output.problem == -1
            % This is the dreaded unknown state from mosek. Try without
            % objective to see if it is infeasible?
            ptest = relaxed_p;
            ptest.c = ptest.c*0;ptest.Q = ptest.Q*0;
            outputtest = bnb_solvelower(lowersolver,ptest,upper,lower,x_min,allSolutions);
            if outputtest.problem == 1
                output.problem = 1;
            else
                output.problem = -1;
            end
        end                       
        
        if output.problem == 9
            unknownErrorCount = unknownErrorCount + 1;
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
                    [upper1,x_min1] = feval(upperSolversList{k},p,upper,x,poriginal,output,lower);
                    if upper1 < upper
                        x_min = x_min1;
                        allSolutions = [allSolutions x_min1];
                        upper = upper1;
                        if length(stack.nodes)>0
                            [stack,stacklower] = prune(stack,upper,p.options,solved_nodes,p,allSolutions);
                            lower = min(lower,stacklower);
                        end
                        [p,stack] = simpleConvexDiagonalQuadraticPropagation(p,upper,stack);
                    elseif ~isinf(upper1) && upper1 == upper && norm(x_min-x_min1) > 1e-4
                        % Yet another solution with same value
                        allSolutions = [allSolutions x_min1];
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
        if (cost<upper) & feasible
            x_min = x;
            upper = cost;
            allSolutions = [allSolutions x_min];
            [stack,lower] = prune(stack,upper,p.options,solved_nodes,p,allSolutions);            
            [p,stack] = simpleConvexDiagonalQuadraticPropagation(p,upper,stack);
        end
        p = adaptivestrategy(p,upper,solved_nodes);
        keep_digging = 0;
    end
    
    % **************************************
    % Stop digging if it won't give sufficient improvement anyway
    % **************************************
    if cost>upper*(1-p.options.bnb.gaptol)
        keep_digging = 0;
    end
    
    feasibilityHistory(end+1) = feasible;
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
        p0_feasible = 1;
        p1_feasible = 1;
        
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
                
        node1 = newNode(p1,globalindex,'up',TotalIntegerInfeas,TotalBinaryInfeas,1-(x(globalindex)-floor(x(globalindex))),pid);
        pid = pid + 1;
        node0 = newNode(p0,globalindex,'down',TotalIntegerInfeas,TotalBinaryInfeas,1-(x(globalindex)-floor(x(globalindex))),pid);
        pid = pid + 1;
%          
%         for j = 1:length(p.binaryProduct)
%             xx = p.binaryProduct{j}.x;
%             yy = p.binaryProduct{j}.y;
%             if any(p0.ub(xx)==0)
%                 if p0.lb(yy)==1
%                     p0_feasible = 0;
%                     break
%                 else
%                     p0.ub(yy) = 0;
%                 end
%             end
%             if any(p1.ub(xx)==0)
%                  if p1.lb(yy)==1
%                     p1_feasible = 0;
%                     break
%                  else
%                      p1.ub(yy) = 0;
%                  end
%             end
%             if all(p0.lb(xx)==1)
%                  if p0.ub(yy)==0
%                     p0_feasible = 0;
%                     break
%                 else
%                 p0.lb(yy)=1;
%                  end
%             end
%             if all(p1.lb(xx)==1)
%                 if p1.ub(yy)==0
%                     p1_feasible = 0;
%                     break
%                 else
%                 p1.lb(yy)=1;
%                 end
%             end
%             if p0.lb(yy)==1
%                 if any(p0.ub(xx)==0)
%                    p0_feasible = 0;
%                     break  
%                 else
%                   p0.lb(xx) = 1;
%                 end
%             end
%             if p1.lb(yy)==1
%                 if any(p1.ub(xx)==0)
%                    p1_feasible = 0;
%                     break  
%                 else
%                   p1.lb(xx) = 1;
%                 end
%             end
%         end
        
        if ismember(globalindex,p.atmost.variables)
            for j = 1:length(p.atmost.groups)
                xy = p.atmost.groups{j};
                if p.atmost.bounds(j)==1 && any(xy == globalindex)
                    if ~(node0.lb(globalindex)==0 && node0.ub(globalindex)==0)
                        % The variable has been fixed to a non-zero value
                        % Hence, its sister has to be set to 0
                        sisters = xy(xy~=globalindex);
                        for k = sisters
                            if node0.lb(k) > 0 || node0.ub(k) < 0
                                p0_feasible = 0;
                                break
                            else
                                node0.lb(k) = 0;
                                node0.ub(k) = 0;
                            end                            
                        end
                    end
                    if ~(node1.lb(globalindex)==0 && node1.ub(globalindex)==0)
                        sisters = xy(xy~=globalindex);
                        for k = sisters
                            if node1.lb(k) > 0 || node1.ub(k) < 0
                                p1_feasible = 0;
                                break
                            else
                                node1.lb(k) = 0;
                                node1.ub(k) = 0;
                            end                            
                        end
                    end
                end
            end
        end
        
        % Make sure we don't push trivially poor stuff to stack, so reuse
        % pruning code by creating temporary stacks first
        tempstack = stackCreate;
        if p0_feasible            
            tempstack = push(tempstack,node0);
        end
        if p1_feasible
            tempstack = push(tempstack,node1);
        end
        tempstack = prune(tempstack,upper,p.options,solved_nodes,p,allSolutions);
        stack = mergeStack(stack,tempstack);    
    end
        
    if stackLength(stack)>0
        lower = stackLower(stack);
        if can_use_ceil_lower
            lower = ceil(lower-1e-6);
        end
    end
    
    % Close current and proceed to next if gap large
    gap = abs((upper-lower)/(1e-3+abs(upper)+abs(lower)));
    [node,stack] = pull(stack,p.options.bnb.method,x_min,upper);
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
    
    lowerhist = [lowerhist lower];
    upperhist = [upperhist upper];
    stacksizehist = [stacksizehist stackLength(stack)];
      
    if p.options.bnb.verbose;
        if mod(solved_nodes-1,p.options.print_interval)==0 || isempty(node) || (gap == 0) || (lastUpper-1e-6 > upper)
            if p.options.bnb.plot
                hold off
                subplot(1,3,1);               
                l = plot([lowerhist' upperhist']);set(l,'linewidth',2);
                title('Upper/lower bounds')
                subplot(1,3,2);                
                l = plot(stacksizehist);set(l,'linewidth',2);
                title('Open nodes')
                drawnow
                subplot(1,3,3);   
                hist(getStackLowers(stack),25);
                title('Histogram lower bounds')
                drawnow
            end
            if lastUpper > upper
                fprintf(' %4.0f : %12.3E  %7.2f   %12.3E  %2.0f  %8.1f    %s \n',solved_nodes,upper,100*gap,lower,stackLength(stack),etime(clock,bnbsolvertime),'-> Found improved solution!');
            else
                fprintf(' %4.0f : %12.3E  %7.2f   %12.3E  %2.0f  %8.1f    %s \n',solved_nodes,upper,100*gap,lower,stackLength(stack),etime(clock,bnbsolvertime),yalmiperror(output.problem,'',1));
            end
        end
    end
    lastUpper = upper;    
end
if p.options.bnb.verbose;showprogress([num2str2(solved_nodes,3)  ' Finishing.  Cost: ' num2str(upper) ],p.options.bnb.verbose);end
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
        indic = (1 + min(10,abs(p.c(all_variables)))).*(abs(x(all_variables)-round(x(all_variables))));
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
        % index = p.semicont_variables(index);
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

p1.binary_variables = new_binary;%p0.binary_variables;%setdiff1D(p1.binary_variables,variable);
%p1.binary_variables = setdiff(p1.binary_variables,friends);
p1.lower = lower;
p1.depth = p.depth+1;

% % *****************************
% % PROCESS MOST PROMISING FIRST
% % (p0 in top of stack)
% % *****************************
if x(variable)>0.5
    pt=p1;
    p1=p0;
    p0=pt;
end

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

% *****************************
% PROCESS MOST PROMISING FIRST
% *****************************
if lb-current<0.5
    pt=p1;
    p1=p0;
    p0=pt;
end

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

function s = num2str2(x,d,c);
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
    if length(p.integer_variables) == length(p.c) && all(p.c == fix(p.c)) && nnz(p.Q)==0 && isempty(p.evalMap) && nnz(p.variabletype)==0       
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
p.atmost = node.atmost;

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

function [p,stack] = pull(stack,method,x_min,upper)

if stackLength(stack) > 0
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
            j = max(find(depths == i));
            p = getStackNode(stack,j);
            stack = removeStackNode(stack,j);
                     
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
            % Silly for backward testing compatibility. Stack order has
            % changed, to be able to compare some examples, make sure we
            % traverse the tree in exactly the same was as before on ties
            j = max(find(lowers == i));
            p = getStackNode(stack,j);
            stack = removeStackNode(stack,j);
            
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
node1.atmost = p1.atmost;

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
