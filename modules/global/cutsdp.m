function output = cutsdp(p)
% CUTSDP
%
% See also SOLVESDP, BNB, BINVAR, INTVAR, BINARY, INTEGER, LMI

% Author Johan Löfberg
% $Id: cutsdp.m,v 1.10 2010-03-09 14:51:21 joloef Exp $

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
%% ADD CONSTRAINTS 0<x<1 FOR BINARY
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

if p.options.bnb.verbose;            disp(' Node       Infeasibility.     Lower    LP cuts');end;

%% Initialize diagnostic
infeasibility = -inf;
solved_nodes = 0;
feasible = 1;
lower = -inf;
saveduals = 1;

% *************************************************************************
%% Add diagonal cuts to begin with
% *************************************************************************
savedCuts = [];
savedIndicies = [];
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
if p.K.s(1)>0
    top = p.K.f+p.K.l+p.K.q+1;
    for i = 1:length(p.K.s)
        n = p.K.s(i);
        newF=[];
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
                newF = [newF;p.F_struc(top+index-1,:)];
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
while goon
    
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

    % Add lower bound
    if ~isinf(lower)
        p_lp.F_struc = [p_lp.F_struc;-lower p_lp.c'];
        p_lp.K.l = p_lp.K.l + 1;
    end

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

    output = feval(cutsolver,p_lp);
  
    % Remove lower bound (avoid accumulating them)
    if ~isinf(lower)
        p_lp.K.l = p_lp.K.l - 1;
        p_lp.F_struc = p_lp.F_struc(1:end-1,:);
    end
         
    if output.problem == 1 | output.problem == 12
        % LP relaxation was infeasible, hence problem is infeasible
        feasible = 0;
        lower = inf;
        goon = 0;
        x = zeros(length(p.c),1);
        lower = inf;
    else
        % Relaxed solution
        x = output.Primal;
       % load dummy
       % nons = [nons;x'];
       % save dummy nons
        lower = p.f+p.c'*x+x'*p.Q*x;
                        
        infeasibility = 0;
        [p_lp,infeasibility] = add_socp_cut(p,p_lp,x,infeasibility);
        [p_lp,infeasibility] = add_sdp_cut(p,p_lp,x,infeasibility);
                
        if ~isempty(pool)
            res = pool*[1;x];
            j = find(res<0)
            if ~isempty(j)
            p_lp.F_struc = [p_lp.F_struc;pool(j,:)];
            p_lp.K.l = p_lp.K.l + length(j);
            pool(j,:)=[];
            end
        end
            
 
        goon = infeasibility <= p.options.cutsdp.feastol;
        goon = goon & feasible;
        goon = goon & (solved_nodes < p.options.cutsdp.maxiter-1);
    end

    solved_nodes = solved_nodes + 1;
    if p.options.cutsdp.verbose;fprintf(' %4.0f :      %12.3E      %12.3E      %2.0f\n',solved_nodes,infeasibility,lower,p_lp.K.l-p.K.l);end
end

D_struc = [];



function [p_lp,infeasibility] = add_sdp_cut(p,p_lp,x,infeasibility);

if p.K.s(1)>0
    % Add cuts
    top = p.K.f+p.K.l+p.K.q+1;
    for i = 1:1:length(p.K.s)
        n = p.K.s(i);
        X = p.F_struc(top:top+n^2-1,:)*[1;x];
        X = reshape(X,n,n);X = (X+X')/2;
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
                if v(m,m)<-1e-12
                    %d(:,m) = round(d(:,m));
                    %d(:,m) = quadprog(X,zeros(length(X),1),[],[],[],[],-ones(length(X),1),ones(length(X),1),d(:,m));
                    if round(d(:,m))'*X*round(d(:,m)) < 0 %& round(d(:,m))'*X*round(d(:,m))<0.2*d(:,m)'*X*d(:,m)
                         d(:,m) = round(d(:,m));
                       % d(:,m) = quadprog(X,zeros(14,1),[],[],[],[],-ones(14,1),ones(14,1),d(:,m));
                    end
                    bA =  d(:,m)'*(kron(d(:,m),speye(n)).'*p.F_struc(top:top+n^2-1,:));
                    b = bA(:,1);
                    A = -bA(:,2:end);                   
                    newF = real([newF;[b -A]]);                   
                    newcuts = newcuts + 1;
                end
            end
        end
        
        newF(abs(newF)<1e-12) = 0;
        keep=find(any(newF(:,2:end),2));
        newF = newF(keep,:);
        if size(newF,1)>0
            p_lp.F_struc = [p_lp.F_struc;newF];
            p_lp.K.l = p_lp.K.l + size(newF,1);
            [i,j] = sort(p_lp.F_struc*[1;x]);
        end
        top = top+n^2;
    end
end

function [p_lp,infeasibility] = add_socp_cut(p,p_lp,x,infeasibility);

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
        keep=find(any(newF(:,2:end),2));
        newF = newF(keep,:);
        if size(newF,1)>0
            p_lp.F_struc = [p_lp.F_struc;newF];
            p_lp.K.l = p_lp.K.l + size(newF,1);
            [i,j] = sort(p_lp.F_struc*[1;x]);
        end
        top = top+n;
    end
end

function plotyalmippolytope(p_lp)

x = sdpvar(length(p_lp.c),1);
  A = -p_lp.F_struc(p_lp.K.f + (1:p_lp.K.l),2:end);
        b = p_lp.F_struc(p_lp.K.f + (1:p_lp.K.l),1);
        plot([A*x < b,-10<x<10]);