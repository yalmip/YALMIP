function output = cutsdp(p)
% CUTSDP
%
% See also OPTIMIZE, BNB, BINVAR, INTVAR, BINARY, INTEGER

% *************************************************************************
%% INITIALIZE DIAGNOSTICS IN YALMIP
% *************************************************************************
bnbsolvertime = clock;
showprogress('Cutting plane solver started',p.options.showprogress);

% *************************************************************************
%% If we want duals, we may not extract bounds. However, bounds must be
% extracted in discrete problems.
% *************************************************************************
if p.options.cutsdp.recoverdual
    warning('Dual recovery not implemented yet in CUTSDP')
end

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

p.ub = min(p.ub,p.options.cutsdp.variablebound');
p.lb = max(p.lb,-p.options.cutsdp.variablebound');

% *************************************************************************
%% PRE-SOLVE (nothing fancy coded)
% *************************************************************************
if isempty(find(isinf([p.ub;p.lb]))) & p.K.l>0
    [p.lb,p.ub] = tightenbounds(-p.F_struc(1+p.K.f:p.K.f+p.K.l,2:end),p.F_struc(1+p.K.f:p.K.f+p.K.l,1),p.lb,p.ub,p.integer_variables);
end

% *************************************************************************
%% PERTURBATION OF LINEAR COST
% *************************************************************************
p.corig = p.c;
if nnz(p.Q)~=0
    g = randn('seed');
    randn('state',1253); %For my testing, I keep this the same...
    % This perturbation has to be better. Crucial for many real LP problems
    p.c = (p.c).*(1+randn(length(p.c),1)*1e-4);
    randn('seed',g);
end

% *************************************************************************
%% We don't need this
% *************************************************************************
p.options.savesolverinput  = 0;
p.options.savesolveroutput = 0;

% *************************************************************************
%% Display logics
% 0 : Silent
% 1 : Display cut progress
% 2 : Display node solver prints
% *************************************************************************
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
[x_min,solved_nodes,lower,feasible,D_struc] = cutting(p);
%% --

% *************************************************************************
%% CREATE SOLUTION
% *************************************************************************
output.problem = 0;
if ~feasible
    output.problem = 1;
end
if solved_nodes == p.options.cutsdp.maxiter
    output.problem = 3;
end
output.solved_nodes = solved_nodes;
output.Primal       = x_min;
output.Dual = D_struc;
output.Slack = [];
output.solverinput  = 0;
output.solveroutput =[];
output.solvertime   = etime(clock,bnbsolvertime);
%% --

function [x,solved_nodes,lower,feasible,D_struc] = cutting(p)

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

% *************************************************************************
%% Create function handle to solver
% *************************************************************************
cutsolver = p.solver.lower.call;

% *************************************************************************
%% Create copy of model without
%  the SDP part
% *************************************************************************
p_lp = p;
p_lp.F_struc = p_lp.F_struc(1:p.K.l+p.K.f,:);
p_lp.K.s = 0;
p_lp.K.q = 0;

% *************************************************************************
%% DISPLAY HEADER
% *************************************************************************
if p.options.cutsdp.verbose
    disp('* Starting YALMIP cutting plane for MISDP based on MILP');
    disp(['* Lower solver   : ' p.solver.lower.tag]);
    disp(['* Max iterations : ' num2str(p.options.cutsdp.maxiter)]);
end

if p.options.bnb.verbose;
    if p.K.s(1)>0
        disp(' Node       Infeasibility     Lower bound    Upper bound  LP cuts   Infeasible SDP cones');
    else
        disp(' Node       Infeasibility     Lower bound    Upper bound  LP cuts');
    end
end

%% Initialize diagnostic
infeasibility = -inf;
solved_nodes = 0;
feasible = 1;
lower = -inf;
saveduals = 1;


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

goon = 1;
rr = p_lp.integer_variables;
rb = p_lp.binary_variables;
only_solvelp = 0;
pplot = 0;

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

hist_infeasibility = [];
mmm=[];
pool = [];

% Avoid data shuffling later on when creating cuts for SDPs
top = 1+p.K.f + sum(p.K.l)+sum(p.K.q);
% Slicing columns much faster
p.F_struc = p.F_struc';
for i = 1:length(p.K.s)
    p.semidefinite{i}.F_struc = p.F_struc(:,top:top+p.K.s(i)^2-1)';
    p.semidefinite{i}.index = 1:p.K.s(i)^2;
    top = top + p.K.s(i)^2;
end
p.F_struc = p.F_struc';
p.F_struc = p.F_struc(1:p.K.f+p.K.l+sum(p.K.q),:);

upper = inf;
standard_options = p_lp.options;
while goon
    
    p_lp = nodeTight(p,p_lp);
    p_lp = nodeFix(p,p_lp);
    
    % Add lower bound
    if ~isinf(lower)
        p_lp.F_struc = [p_lp.F_struc;-lower p_lp.c'];
        p_lp.K.l = p_lp.K.l + 1;
    end
    
    goon_locally = 1;
    p_lp.options = standard_options;
    while goon_locally
        if p.solver.lower.constraint.inequalities.secondordercone.linear
            ptemp = p_lp;
            ptemp.F_struc = [p_lp.F_struc;p.F_struc(1+p.K.f+p.K.l:p.K.f+p.K.l+sum(p.K.q),:)];
            ptemp.K.q = p.K.q;
            output = feval(cutsolver,ptemp);
        else
            output = feval(cutsolver,p_lp);
        end
        
        % Assume we won't find a feasible solution which we try to improve
        goon_locally = 0;
        % Remove lower bound (avoid accumulating them)
        if ~isinf(lower)
            p_lp.K.l = p_lp.K.l - 1;
            p_lp.F_struc = p_lp.F_struc(1:end-1,:);
        end
        
        infeasible_socp_cones = ones(1,length(p.K.q));
        infeasible_sdp_cones = ones(1,length(p.K.s));
        
        eig_failure = 0;
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
            
            infeasibility = 0;
            [p_lp,infeasibility,infeasible_socp_cones] = add_socp_cut(p,p_lp,x,infeasibility);
            [p_lp,infeasibility,infeasible_sdp_cones,eig_failure] = add_sdp_cut(p,p_lp,x,infeasibility);
            [p_lp,infeasibility] = add_nogood_cut(p,p_lp,x,infeasibility);
            
            if ~isempty(pool)
                res = pool*[1;x];
                j = find(res<0)
                if ~isempty(j)
                    p_lp.F_struc = [p_lp.F_struc;pool(j,:)];
                    p_lp.K.l = p_lp.K.l + length(j);
                    pool(j,:)=[];
                end
            end
            
            if feasible && infeasibility >= p.options.cutsdp.feastol && ~eig_failure
                % This was actually a feasible solution
                upper = min(upper, cost);
                if upper > lower
                    if isa(p_lp.options.cutsdp.resolver,'struct')
                        s = p_lp.options.verbose;
                        p_lp.options = p_lp.options.cutsdp.resolver;
                        p_lp.options.verbose = s;
                        goon_locally = 1;
                    end
                end
            end
            
            goon = infeasibility <= p.options.cutsdp.feastol || output.problem ==3;
            goon = goon & feasible;
            goon = goon || eig_failure;% not psd, but no interesting eigenvalue correctly computed
            goon = goon & (solved_nodes < p.options.cutsdp.maxiter-1);
            goon = goon & ~(upper <=lower);
        end
        
        solved_nodes = solved_nodes + 1;
        if eig_failure
            infeasibility = nan;
        end
        if p.options.cutsdp.verbose;
            if p.K.s(1)>0
                if output.problem == 3
                    fprintf(' %4.0f :    %12.3E      %12.3E*  %12.3E     %2.0f         %2.0f/%2.0f\n',solved_nodes,infeasibility,lower,upper,p_lp.K.l-p.K.l,nnz(infeasible_sdp_cones),length(p.K.s));
                else
                    fprintf(' %4.0f :    %12.3E      %12.3E   %12.3E     %2.0f         %2.0f/%2.0f\n',solved_nodes,infeasibility,lower,upper,p_lp.K.l-p.K.l,nnz(infeasible_sdp_cones),length(p.K.s));
                end
            else
                fprintf(' %4.0f :      %12.3E      %12.3E   %12.3E     %2.0f\n',solved_nodes,infeasibility,lower,upper,p_lp.K.l-p.K.l);
            end
        end
        
    end
    
end

D_struc = [];



function [p_lp,worstinfeasibility,infeasible_sdp_cones,eig_computation_failure] = add_sdp_cut(p,p_lp,x,infeasibility_in);

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
    %    lin = p_lp.K.l;
        while iter <= p.options.cutsdp.maxprojections & (infeasibility(end) < -p.options.cutsdp.feastol) && keep_projecting
            % Add one cut b + a'*x >= 0 (if x infeasible)
            %l0 = p_lp.K.l;
            [X,p_lp,infeasibility(iter),a,b,failure] = add_one_sdp_cut(p,p_lp,x,i);                        
            eig_computation_failure = eig_computation_failure & failure;
            if infeasibility(iter) < p_lp.options.cutsdp.feastol && p.options.cutsdp.cutlimit > 0
                % Project current point on the hyper-plane associated with
                % the most negative eigenvalue and move towards the SDP
                % feasible region, and the iterate a couple of iterations
                % to generate a deeper cut
                x0 = x;
                try
                    x = x + a*(-b-a'*x)/(a'*a);
                catch
                    
                end
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


function  [X,p_lp,infeasibility,asave,bsave,failure] = add_one_sdp_cut(p,p_lp,x,i);

newcuts = 0;
newF = [];
n = p.K.s(i);
X = p.semidefinite{i}.F_struc*sparse([1;x]);
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
    [R,indefinite] = chol(X+eye(length(X))*1e-12);
    if indefinite
        infeasibility = -1;
    else
        infeasibility = 0;
    end
    return
end

% For not too large problems, we simply go with a dense
% eigenvalue/vector computation
if  0%n <= p_lp.options.cutsdp.switchtosparse
    [d,v] = eig(X);
    failure = 0;
else
    % Try to perform a block-diagonalization of the current solution,
    % and compute eigenvalue/vectorsa for each block.
    % Sparse eigenvalues can easily fails so we catch info about this
    [d,v,permutation,failure] = dmpermblockeig(X,p_lp.options.cutsdp.switchtosparse);
end
if ~isempty(v)
    d(abs(d)<1e-12)=0;
    infeasibility = min(diag(v));
else
    infeasibility = 0;
end

if infeasibility<0
    [ii,jj] = sort(diag(v));
    
    if ~isempty(permutation)
        [~,inversepermutation] = ismember(1:length(permutation),permutation);
        localFstruc = p.semidefinite{i}.F_struc';
    else
        localFstruc = p.semidefinite{i}.F_struc';
    end
    
    for m = jj(1:min(length(jj),p.options.cutsdp.cutlimit))'
        if v(m,m)<-1e-12
            if 0
                index = reshape(1:n^2,n,n);
                indexpermuted = index(permutation,permutation);
                indexused = index(find(d(:,m)),find(d(:,m)));
                localFstruc = p.F_struc(indexused,:);
                dd=d(find(d(:,m)),m);
                bA = dd'*(kron(dd,speye(length(dd))).'*localFstruc);
            else
                try
                    if ~isempty(permutation)
                        dhere = sparse(d(inversepermutation,m));
                    else
                        dhere = sparse(d(:,m));
                    end
                    if nnz(dhere)>100
                        [~,ii] = sort(-abs(dhere));
                        dhere(abs(dhere) <= abs(dhere(ii(100))))=0;
                    end
                    dd = dhere*dhere';dd = dd(:);
                    index = p.semidefinite{i}.index;
                    used = find(dd);
                    bA = (localFstruc(:,index(used))*sparse(dd(used)))';
                catch
                    bA = sparse(dd(used))'*p.F_struc(index(used),:);
                end
            end
            b = bA(:,1);
            A = -bA(:,2:end);
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

newF(abs(newF)<1e-12) = 0;
keep=find(any(newF(:,2:end),2));
newF = newF(keep,:);
if size(newF,1)>0
    p_lp.F_struc = [p_lp.F_struc(1:p_lp.K.f,:);newF;p_lp.F_struc(1+p_lp.K.f:end,:)];
    p_lp.K.l = p_lp.K.l + size(newF,1);
end

function [p_lp,infeasibility] = add_nogood_cut(p,p_lp,x,infeasibility)
if length(x) ==  length(p.binary_variables)
    % Add a nogood cut. Might already have been generated by
    % the SDP cuts, but it doesn't hurt to add it
    zv = find(x == 0);
    nz = find(x == 1);
    a = zeros(1,length(x));
    a(zv) = 1;
    a(nz) = -1;
    b = length(x)-length(zv)-1;
    newF = [b a];
    p_lp.F_struc = [p_lp.F_struc(1:p_lp.K.f,:);newF;p_lp.F_struc(1+p_lp.K.f:end,:)];
    p_lp.K.l = p_lp.K.l + 1;
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
plot([A*x <= b, p.lb <= x <= p.ub],x,'b',[],sdpsettings('plot.shade',.2));
