function output = cutsdp(p)
%CUTSDP         Cutting plane solver for mixed integer SDP problems
%
% CUTSDP is never called by the user directly, but is called by
% YALMIP from OPTIMIZE, by choosing the solver tag 'cutsdp' in sdpsettings.
%
% The behaviour of CUTSDP can be altered using the fields
% in the field 'cutsdp' in SDPSETTINGS
%
% solver           Solver for the relaxed problems (standard solver tag, see SDPSETTINGS)
%
% maxiter          Maximum number of nodes explored
%
% maxtime          Maximum time allowed
%
% feastol          Tolerance for declaring constraints as feasible
%
% cutlimit         Maximum number of LP cuts added in every iteration
%
% twophase         Start by solving integer-relaxed LPs for a while
%
% See also OPTIMIZE, BNB, BINVAR, INTVAR, BINARY, INTEGER

% *************************************************************************
%% INITIALIZE DIAGNOSTICS IN YALMIP
% *************************************************************************
bnbsolvertime = clock;
showprogress('Cutting plane solver started',p.options.showprogress);
optionsIn = p.options;

% ********************************
%% Remove options if none has been changed
%% Improves peroformance when calling solver many times 
%% (some solvers are really slow at checking options)
% ********************************
p.options = pruneOptions(p.options);

% ********************************
%% Always try to warm-start
% ********************************
p.options.usex0 = 1;

% Arg, new format
p.options.cutsdp.feastol = -p.options.cutsdp.feastol;

% *************************************************************************
%% If we want duals, we may not extract bounds. However, bounds must be
% extracted in discrete problems.
% *************************************************************************
p.noninteger_variables = setdiff(1:length(p.c),[p.integer_variables p.binary_variables p.semicont_variables]);

% *************************************************************************
%% Define infinite bounds
% *************************************************************************
if isempty(p.ub)
    p.ub = repmat(inf,length(p.c),1);
end
if isempty(p.lb)
    p.lb = repmat(-inf,length(p.c),1);
end

% *************************************************************************
%% ADD CONSTRAINTS 0<=x<=1 FOR BINARY
% *************************************************************************
if ~isempty(p.binary_variables)
    p.ub(p.binary_variables) =  min(p.ub(p.binary_variables),1);
    p.lb(p.binary_variables) =  max(p.lb(p.binary_variables),0);
end

% *************************************************************************
%% Extract better bounds from model
% *************************************************************************
if ~isempty(p.F_struc)
    [lb,ub,used_rows_eq,used_rows_lp] = findulb(p.F_struc,p.K);
    if ~isempty([used_rows_eq(:);used_rows_lp(:)])
        lower_defined = find(~isinf(lb));
        if ~isempty(lower_defined)
            p.lb(lower_defined) = max(p.lb(lower_defined),lb(lower_defined));
        end
        upper_defined = find(~isinf(ub));
        if ~isempty(upper_defined)
            p.ub(upper_defined) = min(p.ub(upper_defined),ub(upper_defined));
        end
        p.F_struc(p.K.f+used_rows_lp,:)=[];
        p.F_struc(used_rows_eq,:)=[];
        p.K.l = p.K.l - length(used_rows_lp);
        p.K.f = p.K.f - length(used_rows_eq);
    end
end

% *************************************************************************
%% ADD CONSTRAINTS 0<x<1 FOR BINARY
% *************************************************************************
if ~isempty(p.binary_variables)
    p.ub(p.binary_variables) =  min(p.ub(p.binary_variables),1);
    p.lb(p.binary_variables) =  max(p.lb(p.binary_variables),0);
end

% *************************************************************************
%% PRE-SOLVE (nothing fancy coded)
% *************************************************************************
if isempty(find(isinf([p.ub;p.lb]))) & p.K.l>0
    [p.lb,p.ub] = tightenbounds(-p.F_struc(1+p.K.f:p.K.f+p.K.l,2:end),p.F_struc(1+p.K.f:p.K.f+p.K.l,1),p.lb,p.ub,p.integer_variables);
end

% *************************************************************************
%% We don't need this
% *************************************************************************
p.options.savesolverinput  = 0;
p.options.savesolveroutput = 0;

% *************************************************************************
%% Display logics
% 0 : Silent
% <1: display every 1/N th iteration
% 1 : Display cut progress
% 3 : Display node solver prints
% *************************************************************************
if p.options.verbose ~= fix(p.options.verbose)
    p.options.print_interval = ceil(1/p.options.verbose);
    p.options.verbose = ceil(p.options.verbose);
else
    p.options.print_interval = 1;
end
switch max(min(p.options.verbose,3),0)
    case 0
        p.options.cutsdp.verbose = 0;
    case 1
        p.options.cutsdp.verbose = 1;
        p.options.verbose = 0;
    case 2
        p.options.cutsdp.verbose = 2;
        p.options.verbose = 0;
    case 3
        p.options.cutsdp.verbose = 2;
        p.options.verbose = 1;
    otherwise
        p.options.cutsdp.verbose = 0;
        p.options.verbose = 0;
end

% *************************************************************************
%% START CUTTING
% *************************************************************************
cutsdpsolvertime = clock;
[x_min,solved_nodes,lower,feasible,D_struc,interrupted] = cutting(p);

% *************************************************************************
%% CREATE SOLUTION
% *************************************************************************
if interrupted
    output.problem = 16;
else
    output.problem = 0;
    if ~feasible
        output.problem = 1;
    end
    if solved_nodes == p.options.cutsdp.maxiter
        output.problem = 3;
    elseif etime(clock,cutsdpsolvertime) > p.options.cutsdp.maxtime
        output.problem = 3;
    end
end
output.solved_nodes = solved_nodes;
output.Primal       = x_min;
output.Dual = D_struc;
output.Slack = [];
output.solverinput  = 0;
if optionsIn.savesolveroutput
    output.solveroutput.solved_nodes =solved_nodes;
else
    output.solveroutput =[];
end
output.solvertime   = etime(clock,bnbsolvertime);
%% --

function [x,solved_nodes,lower,feasible,D_struc,interrupted] = cutting(p)

interrupted = 0;
% *************************************************************************
%% Sanity check
% *************************************************************************
if any(p.lb>p.ub)
    x = zeros(length(p.c),1);
    solved_nodes = 0;
    lower = inf;
    feasible = 0;
    D_struc = [];
    return
end

cutsdpsolvertime = clock;

% Solver feasibility too low can cause issues in the termination. 
% However, some solver installations will fail if we do this, so it is not
% default (we have way to combat low precision anyway)
if p.options.cutsdp.adjustsolvertol
    p = adjustSolverPrecision(p);
end

% *************************************************************************
%% Create function handle to solver
% *************************************************************************
cutsolver = p.solver.lower.call;

p = addImpliedSDP(p);

% *************************************************************************
%% Create copy of model without the Conic part
% *************************************************************************
p_lp = p;
p_lp.F_struc = p_lp.F_struc(1:p.K.l+p.K.f,:);
p_lp.K.s = 0;
p_lp.K.q = 0;

p_original = p;

% *************************************************************************
%% DISPLAY HEADER
% *************************************************************************
if p.options.cutsdp.verbose
    disp('* Starting YALMIP cutting plane for MISDP based on MILP');
    disp(['* Lower solver   : ' p.solver.lower.tag]);
    disp(['* Max iterations : ' num2str(p.options.cutsdp.maxiter)]);
    disp(['* Max time       : ' num2str(p.options.cutsdp.maxtime)]);
end

if p.options.cutsdp.verbose
    if p.K.s(1)>0
        disp(' Node    Phase        Cone infeas   Integrality infeas  Lower bound  Upper bound  LP cuts  Elapsed time');
    else
        disp(' Node    Phase        Cone infeas   Integrality infeas  Lower bound  Upper bound  LP cuts  Elapsed time');
    end
end

% Rhs of SOCP has to be non-negative
if ~p.solver.lower.constraint.inequalities.secondordercone.linear
    p_lp = addSOCPCut(p,p_lp);
end

% SDP diagonal has to be non-negative
p_lp = addDiagonalCuts(p,p_lp);

% Experimentation with activation cuts on 2x2 structures in problems with
% all binary variables 2x2 = constant not psd + M(x) means some x has to be
% non-zero
p_lp = addActivationCuts(p,p_lp);
p_lp = removeRedundant(p_lp);

only_solvelp = 0;

% *************************************************************************
% Crude lower bound
% FIX for quadratic case
% *************************************************************************
lower = 0;
if nnz(p.Q) == 0
    for i = 1:length(p.c)
        if p.c(i)>0
            if isinf(p.lb(i))
                lower = -inf;
                break
            else
                lower = lower + p.c(i)*p.lb(i);
            end
        elseif p.c(i)<0
            if isinf(p.ub(i))
                lower = -inf;
                break
            else
                lower = lower + p.c(i)*p.ub(i);
            end
        end
    end
end
%lower = sum(sign(p.c).*(p.lb));
if isinf(lower) | nnz(p.Q)~=0
    lower = -1e6;
end

% *************************************************************************
% Experimental stuff for variable fixing
% *************************************************************************
if p.options.cutsdp.nodefix & (p.K.s(1)>0)
    top=1+p.K.f+p.K.l+sum(p.K.q);
    for i=1:length(p.K.s)
        n=p.K.s(i);
        for j=1:size(p.F_struc,2)-1;
            X=full(reshape(p.F_struc(top:top+n^2-1,j+1),p.K.s(i),p.K.s(i)));
            X=(X+X')/2;
            v=real(eig(X+sqrt(eps)*eye(length(X))));
            if all(v>=0)
                sdpmonotinicity(i,j)=-1;
            elseif all(v<=0)
                sdpmonotinicity(i,j)=1;
            else
                sdpmonotinicity(i,j)=nan;
            end
        end
        top=top+n^2;
    end
else
    sdpmonotinicity=[];
end

% Avoid data shuffling later on when creating cuts for SDPs
top = 1 + p.K.l+p.K.f+sum(p.K.q);
% Slicing columns much faster
p.F_struc = p.F_struc';
for i = 1:length(p.K.s)
    p.semidefinite{i}.F_struc = p.F_struc(:,top:top+p.K.s(i)^2-1)';
    if nnz(p.semidefinite{i}.F_struc)/numel(p.semidefinite{i}.F_struc) > .2
        p.semidefinite{i}.F_struc = full(p.semidefinite{i}.F_struc);
    end
    p.semidefinite{i}.index = 1:p.K.s(i)^2;
    top = top + p.K.s(i)^2;
end
p.F_struc = p.F_struc';
p.F_struc = p.F_struc(1:p.K.f+p.K.l+sum(p.K.q),:);
p.sdpsymmetry = [];
p = detect3x3SymmetryGroups(p);
[p,p_lp] = addSymmetryCuts(p,p_lp);

integer_variables = [p.binary_variables p.integer_variables];

standard_options = p_lp.options;
if p.options.cutsdp.twophase || isempty(integer_variables)
    integerPhase = 0;
else
    integerPhase = 1;
end

% Some structures for keeping a pool of cuts
% silly but we include equalities for simplicity, always active of course)
activity = zeros(p_lp.K.f + p_lp.K.l,1);
prevRem = [];

%% Initialize
x = [];
saveduals = 1;
infeasibility = -inf;
solved_nodes = 0;
feasible = 1;
lower = -inf;
upper = inf;
lowerHistory = -inf;
phaseHistory = [];
goon = 1;
x_found_by_sdp_pump = [];
p_original.sdpPumpData = [];
p_lp_unused = p_lp;
p_lp_unused.F_struc = [];
p_lp_unused.K.f = 0;
p_lp_unused.K.l = 0;
while goon

    % Keep history of what we have been doing. Used for some diagnostics
    phaseHistory = [phaseHistory integerPhase ];
       
    % Basic presolve etc (currently not used)
    p_lp = nodeTight(p,p_lp);
    p_lp = nodeFix(p,p_lp);
        
    % Add upper bound constraint, speeds up the MILP solver slightly somtimes    
    if ~isinf(upper) && (nnz(p_lp.Q)==0)
        p_lp = addLinearCut(p_lp,[upper -p_lp.c']);     
    end
    
    p_lp.options = standard_options;
    
    % Solve relaxed problem
    ptemp = p_lp;
    ptemp.x0 = x;
    if p.solver.lower.constraint.inequalities.secondordercone.linear
        ptemp = p_lp;
        ptemp.F_struc = [p_lp.F_struc;p.F_struc(1+p.K.f+p.K.l:p.K.f+p.K.l+sum(p.K.q),:)];
        ptemp.K.q = p.K.q;
    end
    ptemp = adjustMaxTime(ptemp,ptemp.options.cutsdp.maxtime,etime(clock,cutsdpsolvertime));
    if integerPhase
        output = feval(cutsolver,ptemp);
        if output.problem == 12
            ptemp.c = ptemp.c*0;
            ptemp.Q = ptemp.Q*0;
            output = feval(cutsolver,ptemp);
            if output.problem == 0
                output.problem = 2;
            end
        end
        if min(ptemp.F_struc*[1;output.Primal]) < -abs(p.options.cutsdp.feastol)
            % Ugly hack
            ptemp.F_struc = 1000*ptemp.F_struc;
            output = feval(cutsolver,ptemp);
            ptemp.F_struc = 0.001*ptemp.F_struc;
        end
    else
        ptemp.binary_variables = [];
        ptemp.integer_variables = [];
        output = feval(cutsolver,ptemp);
        if output.problem == 12
            ptemp.c = ptemp.c*0;
            ptemp.Q = ptemp.Q*0;
            output = feval(cutsolver,ptemp);
            if output.problem == 0
                output.problem = 2;
            end
        end
    end
     
    % Remove upper bounds if we added those (avoid accumulating them)   
    if ~isinf(upper) && (nnz(p_lp.Q)==0)
        p_lp.K.l = p_lp.K.l - 1;
        p_lp.F_struc = p_lp.F_struc(1:end-1,:);
    end
       
    % Detect inactive cuts during integer-relaxed phase. Drop these, but
    % same them in a pool to add them later if required
    % [p_lp,activity,prevRem] = handleCutPool(p_lp,activity,prevRem,output,integerPhase,15);
    
    infeasible_socp_cones = ones(1,length(p.K.q));
    infeasible_sdp_cones = ones(1,length(p.K.s));
    
    eig_failure = 0;
    currentPhase = integerPhase ;
    if output.problem == 1 | output.problem == 12
        % LP relaxation was infeasible, hence problem is infeasible
        feasible = 0;
        lower = inf;
        goon = 0;
        x = zeros(length(p.c),1);
        lower = inf;
        cost = inf;
    else
        % Relaxed solution
        x = output.Primal;
        cost = p.f+p.c'*x+x'*p.Q*x;
        if output.problem == 0
            lower = cost;
        end
        
        was_lp_really_feasible = checkfeasiblefast(ptemp,x,-p.options.cutsdp.feastol);               
        
       if p.options.cutsdp.sdppump & integerPhase                     
            [xtemp,upptemp,pumpPossible,p_original,pumpSuccess] = sdpPump(p_original,x,-p.options.cutsdp.feastol);                                               
        else
            upptemp = inf;
        end
        
        infeasibility = 0;
        [p_lp,infeasibility,infeasible_socp_cones] = add_socp_cut(p,p_lp,x,infeasibility);        
        [p_lp,infeasibility,infeasible_sdp_cones,eig_failure] = add_sdp_cut(p,p_lp,x,infeasibility,p_original); 
                        
        % Don't add no-good if we time out etc, for purely integer problems
        % or problems where we analytically can compute the continuous from
        % given integers
        if was_lp_really_feasible && isinf(upptemp) && (output.problem == 0 && ((integerPhase && pumpPossible) || (integerPhase && (length(p.integer_variables)==length(p.c)) && output.problem == 0))            )
            p_lp = add_nogood_cut(p,p_lp,x,infeasibility);                    
        end
        
        if output.problem == 0 && ((integerPhase && pumpPossible) || (integerPhase && (length(p.integer_variables)==length(p.c)) && output.problem == 0))
          p_lp = add3x3sdpsymmetrycut(p,p_lp,x);                   
        end
        
        if upptemp < upper
             upper = upptemp;
             x_found_by_sdp_pump = xtemp;
         end
                    
        if p.options.cutsdp.plot
            plotP(p_lp);drawnow
        end
        
        % Track improvement in lower bound. We terminate the
        % integrality-relaxed phase if no objective improvement
        lowerHistory = [lowerHistory lower];
        if length(lowerHistory) > 10
            b =  mean(lowerHistory(end-10:end-1))+0.01*abs(mean(lowerHistory(end-10:end-1)))            ;
            noprogress = lower <= b;
        else
            noprogress = 0;
        end
               
        sdpInfeasibility = infeasibility;
        % We terminate integrality-relaxed phase if solution is close to
        % SDP-feasible, or no progress in objective, or too many nodes
        if ~isempty(integer_variables) && (infeasibility >= p.options.cutsdp.feastol || solved_nodes > 2000 || noprogress) && ~integerPhase
            integerPhase = 1;
            infeasibility = infeasibility  - inf;
            feasible = 1;          
        elseif infeasibility >= p.options.cutsdp.feastol && integerPhase
            feasible = 1;
            upper = cost;
        end
        
        goon = infeasibility <= p.options.cutsdp.feastol || output.problem ==3 || output.problem == 2;
        
        if length(lowerHistory)>1 && phaseHistory(end)==0
            % We made a massive jump when turning on integrality. This
            % solution is not affected by the possible numerical problems
            % asscociated with adding a lower bound to the relaxed problem
            bigImprovementbyInteger = phaseHistory(end)==1 && phaseHistory(end-1)==0 && lowerHistory(end) >= lowerHistory(end-1) + 0.01*(1e-3 + lowerHistory(end-1));
        elseif length(lowerHistory)==1
            bigImprovementbyInteger = 1;
        else
            bigImprovementbyInteger = 0;
        end
                              
        goon = goon & feasible;
        goon = goon || eig_failure;% not psd, but no interesting eigenvalue correctly computed
        goon = goon & (solved_nodes < p.options.cutsdp.maxiter-1);
        goon = goon & ~(upper <=lower);      
        goon = goon && lower < upper;                       
        if ~isinf(upper) && ~isinf(lower)
            gap = abs((upper-lower)/(1e-3+abs(upper)+abs(lower)));
        else
            gap = inf;
        end
        goon = goon && gap >= p.options.cutsdp.gaptol;        
        goon = goon && (etime(clock,cutsdpsolvertime) < p.options.cutsdp.maxtime);
        goon = goon && ~(output.problem == 16);
    end
    
    solved_nodes = solved_nodes + 1;
    if eig_failure
        infeasibility = nan;
    end
    
    phaseString = {'Continuous','Integer   '};
    integerInfeasibility = sum(abs(x(p.integer_variables)-round(x(p.integer_variables)))) + sum(abs(x(p.binary_variables)-round(x(p.binary_variables))));
    if p.options.cutsdp.verbose
        if mod(solved_nodes-1,p.options.print_interval)==0 || goon == 0
            fprintf(' %4.0f :  %s  %11.3E %12.3E      %14.3E   %11.3E     %3.0f     %5.1f\n',solved_nodes,phaseString{currentPhase+1},sdpInfeasibility,integerInfeasibility,lower,upper,p_lp.K.l-p.K.l,etime(clock,cutsdpsolvertime));
        end
    end
end
% We perhaps terminated before the LP problems found a feasible solution,
% but the SDP pump found a solution along the way, so return that feasible
% solution.
if ~checkfeasiblefast(p_lp,x,-p.options.cutsdp.feastol)
    if ~isempty(x_found_by_sdp_pump)
        x = x_found_by_sdp_pump;
    end
end
D_struc = [];
if output.problem == 16
    interrupted = 1;
end


function [p_lp,worstinfeasibility,infeasible_sdp_cones,eig_computation_failure] = add_sdp_cut(p,p_lp,x,infeasibility_in,p_original);

worstinfeasibility = infeasibility_in;
eig_computation_failure = 0;
infeasible_sdp_cones = zeros(1,length(p.K.s));
if p.K.s(1)>0
    % Solution found by MILP solver
    xsave = x;
    infeasibility = -1;
    eig_computation_failure = 1;
    for i = 1:1:length(p.K.s)
        x = xsave;
        iter = 1;
        keep_projecting = 1;
        infeasibility = 0;                
        psave = p_lp;
        while iter <= p.options.cutsdp.maxprojections+1 & (infeasibility(end) < -p.options.cutsdp.feastol) && keep_projecting
            % Add cuts b + a'*x >= 0 (if x infeasible)            
            [X,p_lp,infeasibility(iter),a,b,failure] = add_one_sdp_cut(p,p_lp,x,i,p_original);
               
            eig_computation_failure = eig_computation_failure & failure;
            if ~isempty(a) && infeasibility(iter) < p_lp.options.cutsdp.feastol && p.options.cutsdp.cutlimit > 0
                % Project current point on the hyper-plane associated with
                % the most negative eigenvalue and move towards the SDP
                % feasible region, and then iterate a couple of iterations
                % to generate a deeper cut                
                x0 = x;                
                x = x + a*(-b-a'*x)/(a'*a);
                x(x <= p.lb) = p.lb(x <= p.lb);                
                x(x >= p.ub) = p.ub(x >= p.ub);               
                keep_projecting = norm(x-x0)>= p.options.cutsdp.projectionthreshold;
            else
                keep_projecting = 0;
            end
            worstinfeasibility = min(worstinfeasibility,infeasibility(iter));
            iter = iter + 1;
        end
        infeasible_sdp_cones(i) =  infeasibility(1) < p_lp.options.cutsdp.feastol;
    end
else
    worstinfeasibility = min(worstinfeasibility,0);
end


function  [X,p_lp,infeasibility,asave,bsave,failure] = add_one_sdp_cut(p,p_lp,x,i,p_original);

newcuts = 0;
newF = [];
n = p.K.s(i);
if numel(x)/length(x) < .1
    X = p.semidefinite{i}.F_struc*sparse([1;x]);
else
    X = p.semidefinite{i}.F_struc*[1;x];
end
X = reshape(X,n,n);X = (X+X')/2;
asave = [];
bsave = [];
% First check if it happens to be psd. Then we are done. Quicker
% than computing all eigenvalues
% This also acts as a slight safe-guard in case the sparse eigs
% fails to prove that the smallest eigenvalue is non-negative
%[R,indefinite] = chol(X+eye(length(X))*1e-12);
%if indefinite

% User is trying to solve by only generating no-good cuts
permutation = [];
failure = 0;
if p.options.cutsdp.cutlimit == 0
     [d,v] = eig(full(X));
     infeasibility = v(1,1);   
    return
end

% For not too large problems, we simply go with a dense
% eigenvalue/vector computation
if  n <= p_lp.options.cutsdp.switchtosparse
    [d,v] = eig(full(X));
    failure = 0;
else
    % Try to perform a block-diagonalization of the current solution,
    % and compute eigenvalue/vectorsa for each block.
    % Sparse eigenvalues can easily fails so we catch info about this
    [d,v,permutation,failure] = dmpermblockeig(X,p_lp.options.cutsdp.switchtosparse);
end
d(abs(d)<1e-12)=0;
infeasibility = min(diag(v));
if infeasibility<0
    [ii,jj] = sort(diag(v));
    
    if ~isempty(permutation)
        [~,inversepermutation] = ismember(1:length(permutation),permutation);      
    end
    
    for m = jj(1:min(length(jj),p.options.cutsdp.cutlimit))'
        if v(m,m)<=-1e-12                       
            try
                if ~isempty(permutation)
                    dhere = d(inversepermutation,m);
                else
                    dhere = d(:,m);
                end
                dd = dhere*dhere';dd = dd(:);
                bA = dd'*p.semidefinite{i}.F_struc;
                if numel(bA)/nnz(bA) < .1
                    bA = sparse(bA);
                end
            end
            b = bA(:,1);
            A = -bA(:,2:end);            
            if isempty(p_lp.F_struc) || ~any(sum(abs(p_lp.F_struc-repmat([b -A],size(p_lp.F_struc,1),1)),2)<= 1e-12)
                newF = real([newF;[b -A]]);
                newcuts = newcuts + 1;
                if isempty(asave)
                    A(abs(A)<1e-12)=0;
                    b(abs(b)<1e-12)=0;
                    asave = -A(:);
                    bsave = b;
                end
            end
        end
    end
end

newF(abs(newF)<1e-12) = 0;
keep=find(any(newF(:,2:end),2));
newF = newF(keep,:);
if size(newF,1)>0
    newF(:,1) = newF(:,1) + 0*0.02*abs(p_lp.options.cutsdp.feastol);
    p_lp.F_struc = [p_lp.F_struc(1:p_lp.K.f,:);p_lp.F_struc(1+p_lp.K.f:end,:);newF];
    p_lp.K.l = p_lp.K.l + size(newF,1);
end

function [p_lp] = add_nogood_cut(p,p_lp,x,infeasibility)
% Add a nogood cut. Might already have been generated by
% the SDP cuts, but it doesn't hurt to add it
if length(x) ==  length(p.binary_variables)    
    % Standard binary case
    [b,a] = exclusionCut(x,1);    
    p_lp.F_struc = [p_lp.F_struc(1:p_lp.K.f,:);p_lp.F_struc(1+p_lp.K.f:end,:);b a];
    p_lp.K.l = p_lp.K.l + 1;
elseif all(p_lp.ub(p.integer_variables) <= 0) && all(p_lp.lb(p.integer_variables) >= -1) && length(p.binary_variables)==0
    % Not all variables are integer, but we've analytically showed that
    % this solution isn't possible anyway. Also, this is the negated binary
    % case. TODO: Generalize    
    x(p_lp.noninteger_variables) = [];
    [b,atemp] = exclusionCut(x,-1);
    a = zeros(1,length(p.c));
    a(p.integer_variables) = atemp;
    p_lp.F_struc = [p_lp.F_struc(1:p_lp.K.f,:);p_lp.F_struc(1+p_lp.K.f:end,:);b a];
    p_lp.K.l=p_lp.K.l+1;
end


function [p_lp,infeasibility,infeasible_socp_cones] = add_socp_cut(p,p_lp,x,infeasibility);
infeasible_socp_cones = zeros(1,length(p.K.q));
% Only add these cuts if solver doesn't support SOCP cones
if ~p.solver.lower.constraint.inequalities.secondordercone.linear
    if p.K.q(1)>0
        % Add cuts
        top = p.K.f+p.K.l+1;
        for i = 1:1:length(p.K.q)
            n = p.K.q(i);
            X = p.F_struc(top:top+n-1,:)*[1;x];
            X = [X(1) X(2:end)';X(2:end) eye(n-1)*X(1)];
            Y = randn(n,n);
            newcuts = 1;
            newF = zeros(n,size(p.F_struc,2));
            [d,v] = eig(X);
            infeasibility = min(infeasibility,min(diag(v)));
            dummy=[];
            newF = [];
            if infeasibility<0
                [ii,jj] = sort(diag(v));
                for m = jj(1:min(length(jj),p.options.cutsdp.cutlimit))'%find(diag(v<0))%1:1%length(v)
                    if v(m,m)<0
                        v1 = d(1,m);v2 = d(2:end,m);
                        newF = [newF;p.F_struc(top,:) + 2*v1*v2'*p.F_struc(top+1:top+n-1,:)];
                        newcuts = newcuts + 1;
                    end
                end
            end
            newF(abs(newF)<1e-12) = 0;
            keep= any(newF(:,2:end),2);
            newF = newF(keep,:);
            if size(newF,1)>0
                p_lp.F_struc = [p_lp.F_struc;newF];
                p_lp.K.l = p_lp.K.l + size(newF,1);
                [i,j] = sort(p_lp.F_struc*[1;x]);
            end
            top = top+n;
        end
    end
end

function p_lp = addActivationCuts(p,p_lp)
if p.options.cutsdp.activationcut && p.K.s(1) > 0 && length(p.binary_variables) == length(p.c)
    top = p.K.f + p.K.l+sum(p.K.q)+1;
    for k = 1:length(p.K.s)
        F0 = p.F_struc(top:top+p.K.s(k)^2-1,1);
        %  Fij = p.F_struc(top:top+p.K.s(k)^2-1,2:end);
        %  Fij = sum(Fij | Fij,2);
        F0 = reshape(F0,p.K.s(k),p.K.s(k));
        %  Fij = reshape(Fij,p.K.s(k),p.K.s(k));
        %  Fall = F0 | Fij;
        row = 1;
        added = 0; % Avoid adding more than 2*n cuts (we hope for sparse model...)
        while row <= p.K.s(k)-1 && added <= 2*p.K.s(k)
            % if 1
            j = find(F0(row,:));
            if min(eig(F0(j,j)))<0
                [ii,jj] = find(F0(j,j));
                ii = j(ii);
                jj = j(jj);
                index = sub2ind([p.K.s(k),p.K.s(k)],ii,jj);
                p.F_struc(top + index-1,2:end);
                S = p.F_struc(top + index-1,2:end);
                S = S | S;S = sum(S,1);S = S | S;
                % Some of these have to be different from 0
                p_lp.F_struc = [p_lp.F_struc;-1 S];
                p_lp.K.l = p_lp.K.l + 1;
                
            end
            % else
            j = find(F0(row,:));
            j = j(j>row);
            for col = j(:)'
                if F0(row,row)*F0(col,col)-F0(row,col)^2<0
                    index = sub2ind([p.K.s(k),p.K.s(k)],[row row col],[row col col]);
                    p.F_struc(top + index-1,2:end);
                    S = p.F_struc(top + index-1,2:end);
                    S = S | S;S = sum(S,1);S = S | S;
                    % Some of these have to be different from 0
                    p_lp.F_struc = [p_lp.F_struc;-1 S];
                    added = added + 1;
                    p_lp.K.l = p_lp.K.l + 1;
                end
            end
            row = row + 1;
        end
    end
end

function p_lp = addDiagonalCuts(p,p_lp)
if p.K.s(1)>0
    top = p.K.f+p.K.l+sum(p.K.q)+1;
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        newF=[];
        nouse = [];
        for m = 1:p.K.s(i)
            d = eyev(p.K.s(i),m);
            index = (1+(m-1)*(p.K.s(i)+1));
            ab = p.F_struc(top+index-1,:);
            b =  ab(1);
            a = -ab(2:end);
            % a*x <= b
            pos = find(a>0);
            neg = find(a<0);
            if a(pos)*p.ub(pos) + a(neg)*p.lb(neg)>b
                if length(p.binary_variables) == length(p.c)
                    if all(p.F_struc(top+index-1,2:end) == fix(p.F_struc(top+index-1,2:end)))
                        ab(1) = floor(ab(1));
                        if max(a)<=0 % Exclusive or in disguise
                            ab = sign(ab);
                        end
                    end
                end
                newF = [newF;ab];
            else
                nouse = [nouse m];
            end
        end
        % Clean
        newF(abs(newF)<1e-12) = 0;
        keep=find(any(newF(:,2:end),2));
        newF = newF(keep,:);
        
        p_lp.F_struc = [p_lp.F_struc;newF];
        p_lp.K.l = p_lp.K.l + size(newF,1);
        top = top+n^2;
    end
end

function p_lp = addSOCPCut(p,p_lp)
if p.K.q(1) > 0
    top = p.K.f+p.K.l+1;
    for i = 1:length(p.K.q)
        n = p.K.q(i);
        newF = p.F_struc(top,:);
        
        % Clean
        newF(abs(newF)<1e-12) = 0;
        keep=find(any(newF(:,2:end),2));
        newF = newF(keep,:);
        
        p_lp.F_struc = [p_lp.F_struc;newF];
        p_lp.K.l = p_lp.K.l + size(newF,1);
        top = top+n;
    end
end

function p_lp = nodeTight(p,p_lp);
if p.options.cutsdp.nodetight
    % Extract LP part Ax<=b
    A = -p_lp.F_struc(p_lp.K.f + (1:p_lp.K.l),2:end);
    b = p_lp.F_struc(p_lp.K.f + (1:p_lp.K.l),1);
    c = p_lp.c;
    % Tighten bounds and find redundant constraints
    [p_lp.lb,p_lp.ub,redundant,pss] = milppresolve(A,b,p_lp.lb,p_lp.ub,p.integer_variables,p.binary_variables,ones(length(p.lb),1));
    A(redundant,:) = [];
    b(redundant) = [];
    p_lp.F_struc(p_lp.K.f+redundant,:) = [];
    p_lp.K.l = p_lp.K.l-length(redundant);
end

function p_lp = nodeFix(p,p_lp);
if p.options.cutsdp.nodefix
    % Try to find variables to fix w.l.o.g
    [fix_up,fix_down] = presolve_fixvariables(A,b,c,p_lp.lb,p_lp.ub,sdpmonotinicity);
    p_lp.lb(fix_up)   = p_lp.ub(fix_up);
    p_lp.ub(fix_down) = p_lp.lb(fix_down);
    while ~(isempty(fix_up) & isempty(fix_down))
        [p_lp.lb,p_lp.ub,redundant,pss] = milppresolve(A,b,p_lp.lb,p_lp.ub,p.integer_variables,p.binary_variables,ones(length(p.lb),1));
        A(redundant,:) = [];
        b(redundant) = [];
        p_lp.F_struc(p_lp.K.f+redundant,:) = [];
        p_lp.K.l = p_lp.K.l-length(redundant);
        fix_up = [];
        fix_down = [];
        % Try to find variables to fix w.l.o.g
        [fix_up,fix_down] = presolve_fixvariables(A,b,c,p_lp.lb,p_lp.ub,sdpmonotinicity);
        p_lp.lb(fix_up)   = p_lp.ub(fix_up);
        p_lp.ub(fix_down) = p_lp.lb(fix_down);
    end
end

function p_lp = removeRedundant(p_lp);

F = unique(p_lp.F_struc(1+p_lp.K.f:end,:),'rows');
if size(F,1) < p_lp.K.l
    p_lp.F_struc = [p_lp.F_struc(1:p_lp.K.f,:);F];
    p_lp.K.l = size(F,1);
end


function plotP(p)

b = p.F_struc(1+p.K.f:p.K.f+p.K.l,1);
A = -p.F_struc(1+p.K.f:p.K.f+p.K.l,2:end);
x = sdpvar(size(A,2),1);
plot([A*x <= b, p.lb <= x <= p.ub],x(1:2),'b',[],sdpsettings('plot.shade',.2));


function [p_lp,activity,prevRem] = handleCutPool(p_lp,activity,prevRem,output,integerPhase,activitylimit)
if length(activity) < p_lp.K.l + p_lp.K.f
      activity(p_lp.K.l+p_lp.K.f) = 0;
end
if ~integerPhase
    % Save info about active cuts
    % (for how long have they been consequitively active)
    Inactive_Here = (p_lp.F_struc*[1;output.Primal] >= 1e-3);
    activity(find(Inactive_Here)) = activity(find(Inactive_Here)) + 1;
    activity(find(~Inactive_Here)) = 0;
    % Prune inactive cuts that have been inactive for many iterations
    remove = find(activity(:)>activitylimit);
    prevRem = [prevRem;p_lp.F_struc(remove,:)];
    p_lp.F_struc(remove,:) = [];
    p_lp.K.l = p_lp.K.l - length(remove);
    activity(remove)=[];
end
if ~isempty(prevRem)
    % We're keeping a pool of removed cuts, and if they turn active
    % again, add them and define them as very likely to be kept
    Violated_Here = find((prevRem*[1;output.Primal] <= 0));
    if ~isempty(Violated_Here)
        p_lp = addLinearCut(p_lp,prevRem(Violated_Here,:));
%        p_lp.F_struc = [p_lp.F_struc;prevRem(Violated_Here,:)];
%        p_lp.K.l = p_lp.K.l + length(Violated_Here);
        activity = [activity(:);ones(length(Violated_Here),1)*-20];      
        prevRem(Violated_Here,:)=[];
    end
end

function p = addLinearCut(p,row);
p.F_struc = [p.F_struc;row];
p.K.l = p.K.l + size(row,1);


function p = adjustSolverPrecision(p)

switch p.options.cutsdp.solver
    case 'cplex'
        try
            % Case that user specified, catch error if empty
            p.options.cplex.simplex.tolerances.feasibility = min([abs(p.options.cutsdp.feastol/10) p.options.cplex.simplex.tolerances.feasibility]);
        catch
            % Options has been removed, i.e. user has not specified?
            try
                p.options.cplex.simplex.tolerances.feasibility = min([abs(p.options.cutsdp.feastol/10) p.options.default.cplex.simplex.tolerances.feasibility]);
            catch
                % User is sitting on obsolete cplex/matlab
                disp('WARNING: You should update CPLEX. THe options structure does not work as expected');
            end
        end        
            
    otherwise
end


function p = detect3x3SymmetryGroups(p)

good = zeros(1,length(p.c));
good(p.integer_variables) = 1;
good( (p.lb ~= -1) | (p.ub ~=0)) = 0;
if any(good)
    groups = {};
    for j = 1:length(p.K.s)
        n = 3;
        X = spalloc(p.K.s(j),p.K.s(j),p.K.s(j)^2);
        X(1:n,1:n) = 1;
        index0 = find(X);
        index = index0;
        corner = 0;        
        for block = 1:p.K.s(j)-n
            dataBlock = p.semidefinite{j}.F_struc(index,:);
            used = find(any(dataBlock,1));
            dataBlock = dataBlock(:,used);
            if used(1) == 0;dataBlock = [zeros(size(dataBlock,1),1) dataBlock];end
            v = used;v = v(v>1)-1;
            if all(good(v))
                if isempty(groups)
                    groups{1}.dataBlock = dataBlock;
                    groups{1}.variables{1} = used;
                else
                    found = 0;
                    for k = 1:length(groups)
                        if isequal(length(used),length(groups{1}.variables{1}))
                            if isequal(groups{1}.dataBlock,dataBlock)
                                found = 1;
                                groups{1}.variables{end+1} = used;
                            else
                                % s = groups{1}.dataBlock(:)\dataBlock(:);
                                % norm(s*groups{1}.dataBlock-dataBlock,inf)
                            end
                        end
                    end
                    if ~found
                        groups{end+1}.dataBlock = dataBlock;
                        groups{end}.variables{1} = used;
                    end
                end
            end
            index = index + p.K.s(j)+1;
        end
    end
    if length(groups) > 0
        for i = 1:length(groups)
            if length(groups{i}.variables) > 1
                keep(i) = 1;
            else
                keep(i) = 0;
            end
        end
        groups = {groups{find(keep)}};
        if length(groups) > 0
            for i = 1:length(groups)
                for j = 1:length(groups{i}.variables);
                    v = groups{i}.variables{j};
                    v = v(v>1)-1;
                    groups{i}.variables{j} = v;
                end
            end
        end
    end
else
    groups = {};
end
p.sdpsymmetry = groups;
    
    
function p_lp = add3x3sdpsymmetrycut(p,p_lp,x) 

for j = 1:length(p.sdpsymmetry)    
    excludes = [];
    n = sqrt(size(p.sdpsymmetry{j}.dataBlock,1));
    for i = 1:length(p.sdpsymmetry{j}.variables)        
        if min(eig(reshape(p.sdpsymmetry{j}.dataBlock*[1;x(p.sdpsymmetry{j}.variables{i})],n,n))) < -abs(p_lp.options.cutsdp.feastol)
            excludes = [excludes x(p.sdpsymmetry{j}.variables{i})];
        end
    end
    if ~isempty(excludes)
        newF = [];
        infeasible_combinations = unique(excludes','rows')';
        for k = 1:size(infeasible_combinations,2)
            % Local cut for reduced set
            [b,atemp] = exclusionCut(infeasible_combinations(:,k),-1);
            % Add that cut for every variable groups
            for s = 1:length(p.sdpsymmetry{j}.variables)
                a = spalloc(1,length(p_lp.c),1);
                a(p.sdpsymmetry{j}.variables{s}) = atemp;
                if all(sum(abs(p_lp.F_struc - [b a]),2)<=1e-12)
                    newF = [newF;b a];
                end
            end
        end
        p_lp = addLinearCut(p_lp,newF);
    end
end
    
    
function [p,p_lp] = addSymmetryCuts(p,p_lp)

for j = 1:length(p.sdpsymmetry)    
    if length(p.sdpsymmetry{j}.variables{1}) <= 4
        % We can enumerate easily infeasibles
        excludes = [];
        n = sqrt(size(p.sdpsymmetry{j}.dataBlock,1));
        combs = -dec2decbin(0:2^length(p.sdpsymmetry{j}.variables{1})-1,length(p.sdpsymmetry{j}.variables{1}))';
        for i = 1:size(combs,2)
            if min(eig(reshape(p.sdpsymmetry{j}.dataBlock*[1;combs(:,i)],n,n))) < -abs(p_lp.options.cutsdp.feastol)
                excludes = [excludes combs(:,i)];
            end
        end
        if ~isempty(excludes)
            newF = [];
            infeasible_combinations = unique(excludes','rows')';
            for k = 1:size(infeasible_combinations,2)
                % Local cut for reduced set
                [b,atemp] = exclusionCut(infeasible_combinations(:,k),-1);
                % Add that cut for every variable groups
                for s = 1:length(p.sdpsymmetry{j}.variables)
                    a = spalloc(1,length(p_lp.c),1);
                    a(p.sdpsymmetry{j}.variables{s}) = atemp;                    
                    newF = [newF;b a];                  
                end
            end
            p_lp = addLinearCut(p_lp,newF);
            % Delete, won't need these in the future
            p.sdpsymmetry{j}.dataBlock = [];
            p.sdpsymmetry{j}.variables = [];
        end
    end
end
