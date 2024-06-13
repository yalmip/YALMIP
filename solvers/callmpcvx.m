function output = callmpcvx(p)
%BMIBNB          Branch-and-bound scheme for bilinear programs
%
% BMIBNB is never called by the user directly, but is called by 
% YALMIP from SOLVESDP, by choosing the solver tag 'bmibnb' in sdpsettings
%
% The behaviour of BMIBNB can be altered using the fields
% in the field 'bmibnb' in SDPSETTINGS
%
% WARNING: THIS IS EXPERIMENTAL CODE
%
% bmibnb.lowersolver    - Solver for lower bound [standard solver tag ('')]
% bmibnb.uppersolver    - Solver for upper bound [standard solver tag ('')]
% bmibnb.lpsolver       - Solver for LP bound tightening [standard solver tag ('')]
% bmibnb.branchmethod   - Branch strategy ['maxvol' | 'best' ('best')]
% bmibnb.branchrule     - Branch position ['omega' | 'bisect' ('omega')]
% bmibnb.lpreduce       - Improve variable bounds using LP [ real [0,1] (0 means no reduction, 1 means all variables)
% bmibnb.lowrank        - partition variables into two disjoint sets and branch on smallest [ 0|1 (0)]
% bmibnb.target         - Exit if upper found<target [double (-inf)]
% bmibnb.roottight      - Improve variable bounds in root using full problem [ 0|1 (1)]
% bmibnb.vartol         - Cut tree when x_U-x_L < vartol on all branching variables
% bmibnb.relgaptol      - Tolerance on relative objective error (UPPER-LOWER)/(1+|UPPER|) [real (0.01)]
% bmibnb.absgaptol      - Tolerance on objective error (UPPER-LOWER) [real (0.01)]
% bmibnb.pdtol          - A number is declared non-negative if larger than...[ double (-1e-6)]
% bmibnb.eqtol          - A number is declared zero if abs(x) smaller than...[ double (1e-6)]
% bmibnb.maxiter        - Maximum number of solved nodes [int (100)]
% bmibnb.maxtime        - Maximum CPU time (sec.) [int (3600)]

% ********************************
% INITIALIZE DIAGNOSTICS IN YALMIP
% ********************************
bnbsolvertime = clock; 
showprogress('Branch and bound started',p.options.showprogress);

% *******************************
% Display-logics
% *******************************
switch max(min(p.options.verbose,3),0)
    case 0
        p.options.bmibnb.verbose = 0;
    case 1
        p.options.bmibnb.verbose = 1;
        p.options.verbose = 0;
    case 2
        p.options.bmibnb.verbose = 2;
        p.options.verbose = 0;
    case 3
        p.options.bmibnb.verbose = 2;
        p.options.verbose = 1;
    otherwise 
        p.options.bmibnb.verbose = 0;
        p.options.verbose = 0;    
end

% *******************************
% Actual linear variables
% *******************************
p.linears = find(sum(p.monomtable,2)==1);

% *******************************
% PRE-CALCULATE INDICIES
% FOR BILNEAR VARIABLES
% *******************************
nonlinear = find(~(sum(p.monomtable~=0,2)==1 & sum(p.monomtable,2)==1));
nonlins   = [];
for i = 1:length(nonlinear)
    z = find(p.monomtable(nonlinear(i),:));
    if length(z)==1
        nonlins = [nonlins;nonlinear(i) z z];
    else
        nonlins = [nonlins;nonlinear(i) z(1) z(2)];
    end
end
p.nonlins = nonlins;

p.options.saveduals = 0;

% *******************************
% Select branch variables
% *******************************
p.branch_variables = decide_branch_variables(p);

% ********************************
% Extract bounds from model
% ********************************
if isempty(p.ub)
    p.ub = repmat(inf,length(p.c),1);
end
if isempty(p.lb)
    p.lb = repmat(-inf,length(p.c),1);
end
if ~isempty(p.F_struc)
    [lb,ub,used_rows] = find_lp_bounds(p.F_struc,p.K);
    if ~isempty(used_rows)
        lower_defined = find(~isinf(lb));
        if ~isempty(lower_defined)
            p.lb(lower_defined) = max(p.lb(lower_defined),lb(lower_defined));
        end
        upper_defined = find(~isinf(ub));
        if ~isempty(upper_defined)
            p.ub(upper_defined) = min(p.ub(upper_defined),ub(upper_defined));
        end   
        % Remove linear bounds
        used_rows = used_rows(find(~any(p.F_struc(p.K.f+used_rows,1+nonlinear),2)));
        not_used_rows = setdiff(1:p.K.l,used_rows);
        for i = 1:length(p.KCut.l)
            p.KCut.l(i) = find(not_used_rows==p.KCut.l(i) );
        end
        if ~isempty(used_rows)
            p.F_struc(p.K.f+used_rows,:)=[];
            p.K.l = p.K.l - length(used_rows);
        end       
    end
end
p = clean_bounds(p);
p = updatenonlinearbounds(p);
feasible = all(p.lb<=p.ub);

% ********************************
% Remove empty linear rows
% ********************************
if any(p.K.l)
    empty_rows = find(~any(p.F_struc(p.K.f+1:p.K.f+p.K.l,2:end),2));
    if ~isempty(empty_rows)
        if all(p.F_struc(p.K.f+empty_rows,1)>=0)            
            p.F_struc(p.K.f+empty_rows,:)=[];
            p.K.l = p.K.l - length(empty_rows);         
        else
            feasible = 0;
        end
    end
end

% ********************************
% Tighten bounds at root
% ********************************
if p.options.bmibnb.roottight & feasible
    lowersolver = eval(['@' p.solver.lowercall]);
    c = p.c;
    Q = p.Q;
    mt = p.monomtable;
    p.monomtable = eye(length(c));
    i = 1;
    while i<=length(p.linears) & feasible   
        j = p.linears(i);  
        p.c = eyev(length(p.c),j);
        output = feval(lowersolver,p);  
        if (output.problem == 0) & (output.Primal(j)>p.lb(j))
            p.lb(j) = output.Primal(j);  
            p = updateonenonlinearbound(p,j);
            p = clean_bounds(p);
        end
        if output.problem == 1
            feasible = 0;
        else
            % p = updatenonlinearbounds(p,0,1);
            p.c = -eyev(length(p.c),j);
            output = feval(lowersolver,p);
            if (output.problem == 0) & (output.Primal(j) < p.ub(j))
                p.ub(j) = output.Primal(j);
                p = updateonenonlinearbound(p,j);
                p = clean_bounds(p);
            end    
            if output.problem == 1
                feasible = 0;
            end             
            i = i+1;
        end
    end
    p.lb(p.lb<-1e10) = -inf;
    p.ub(p.ub>1e10) = inf;
    p.c = c;
    p.Q = Q;
    p.monomtable = mt;
end

if feasible
    
    % *******************************
    % Bounded domain?
    % *******************************
    involved = unique(p.nonlins(:,2:3));
    if isinf(p.lb(involved)) | isinf(p.ub(involved))    
        error('You have to bound all complicating variables explicitely (I cannot deduce bounds on all variables)')
        output.Primal = [];
        output.problem = -1;
    end
        
    % *******************************
    % We don't need to save node data
    % *******************************
    p.options.savesolverinput  = 0;
    p.options.savesolveroutput = 0;
    
    % *******************************
    % RUN BRANCH & BOUND
    % *******************************
    [x_min,solved_nodes,lower,upper] = branch_and_bound(p);
    
    % **********************************
    % CREATE SOLUTION
    % **********************************
    output.problem = 0;
    if isinf(upper)
        output.problem = 1;
    end
    if isinf(-lower)
        output.problem = 2;
    end
    if solved_nodes == p.options.bnb.maxiter
        output.problem = 3;
    end
else
    output.problem = 1;
    x_min = repmat(nan,length(p.c),1);
    solved_nodes = 0;
end

output.solved_nodes = solved_nodes;
output.Primal      = x_min;
output.Dual        = [];
output.Slack       = [];
output.infostr      = yalmiperror(output.problem,'BNB');
output.solverinput  = 0;
output.solveroutput =[];
output.solvertime   = etime(clock,bnbsolvertime);

function [x_min,solved_nodes,lower,upper] = branch_and_bound(p)

% ***************************************
% LPs ARE USED IN  BOX-REDUCTION 
% (this is essentially a cutting plane pool)
% ***************************************
p.lpcuts = p.F_struc(1+p.K.f:1:p.K.l,:);

% ***************************************
% Create function handles to solvers
% ***************************************
try
    lowersolver = eval(['@' p.solver.lowercall]); % Local LMI solver
    uppersolver = eval(['@' p.solver.uppercall]); % Local BMI solver
    lpsolver    = eval(['@' p.solver.lpcall]);    % LP solver
catch
    disp(' ');
    disp('The internal branch & bound solver requires MATLAB 6')
    disp('I am too lazy too do the changes to make it compatible')
    disp('with MATLAB 5. If you really need it, contact me...');
    disp(' ');    
    error(lasterr);
end

% ************************************************
% GLOBAL PROBLEM DATA
% ************************************************
c       = p.c;
Q       = p.Q;
f       = p.f;
K       = p.K;
p.options.saveduals = 0;
options = p.options;

% ************************************************
% ORIGINAL PROBLEM (used in LOCAL BMI solver)
% ************************************************
p_upper = p;

% ************************************************
% Remove linear cutting planes from problem
% ************************************************
p_upper.F_struc(p_upper.K.f+p_upper.KCut.l,:)=[];
p_upper.K.l = p_upper.K.l - length(p_upper.KCut.l);

% ************************************************
% Remove sdp cutting planes from problem
% ************************************************
if length(p_upper.KCut.s)>0
    starts = p_upper.K.f+p_upper.K.l + [1 1+cumsum((p_upper.K.s).^2)];
    remove_these = [];
    for i = 1:length(p_upper.KCut.s)
        j = p_upper.KCut.s(i);
        remove_these = [remove_these;(starts(j):starts(j+1)-1)'];
    end
    p_upper.F_struc(remove_these,:)=[];
    p_upper.K.s(p_upper.KCut.s) = [];
end

% ************************************************
% INITILIZATION
% ************************************************
p.depth = 0;
p.dpos = 0; 
p.lower = NaN;
upper   = inf;
lower   = NaN;
gap = inf;
x_min   = zeros(length(p.c),1);
stack   = [];
solved_nodes = 0;
solved_lower = 0;
solved_upper = 0;
solved_lp = 0;

if isempty(p.x0)
    p.x0    = zeros(length(p.c),1);
end

x0 = evaluate_nonlinear(p,p.x0);
upper_residual = resids(p,x0);
x0_feasible = all(upper_residual(1:p.K.f)>=-options.bmibnb.eqtol) & all(upper_residual(1+p.K.f:end)>=options.bmibnb.pdtol);
if p.options.warmstart & x0_feasible
    x_min = x0;
    upper = p.f+p.c'*x0+x0'*Q*x0;
end

% ************************************************
% Branch & bound loop
% ************************************************
if options.bmibnb.verbose>0
    fprintf('******************************************************************************************************************\n')
    fprintf('#node          Was'' up                         gap      upper      node    lower  dpth  stk     Memory Vol-red\n')
    fprintf('******************************************************************************************************************\n')
end

doplot = 0;
if doplot
    close all;
    hold on;
end

t_start = cputime;
go_on  = 1;
while go_on
    
    if doplot;ellipplot(diag([200 50]),1,'y',[p.dpos;-p.depth]);drawnow;end;
    % ********************************************
    % ASSUME THAT WE WON'T FATHOME
    % ********************************************
    keep_digging = 1;
    % ********************************************
    % REDUCE BOX
    % ******************************************** 
    if ~options.bmibnb.lpreduce
        % [p.lb,p.ub] = tightenbounds(-p.F_struc(1+p.K.f:p.K.f+p.K.l,2:end),p.F_struc(1+p.K.f:p.K.f+p.K.l,1),p.lb,p.ub,[]);  
        vol_reduction = 1;
        feasible = 1;
    else     
        [p,feasible,vol_reduction] =  boxreduce(p,upper,lower,lpsolver,options);  
    end
    
    % ********************************************
    % SOLVE LOWER AND UPPER
    % ********************************************
    if feasible               
        
        output = solvelower(p,options,lowersolver); 
        
        info_text = '';
        switch output.problem
            case 1
                if doplot;ellipplot(diag([200 25]),1,'r',[p.dpos;-p.depth]);drawnow;end;
                info_text = 'Infeasible node';
                keep_digging = 0;
                cost = inf;
                feasible = 0;
                
            case 2
                cost = -inf;
                
            case {0,3,4}
                
                x = output.Primal; 
                
                cost = f+c'*x+x'*Q*x;
                
                z = evaluate_nonlinear(p,x);
                p = addsdpcut(p,z);        
                
                % Maybe the relaxed solution is feasible
                relaxed_residual = resids(p_upper,z);
                relaxed_feasible = all(relaxed_residual(1:p.K.f)>=-options.bmibnb.eqtol) & all(relaxed_residual(1+p.K.f:end)>=options.bmibnb.pdtol);
                if relaxed_feasible
                    this_upper = f+c'*z+z'*Q*z;
                    if (this_upper < (1-1e-5)*upper) & (this_upper < upper - 1e-5)
                        x_min = z;
                        upper = this_upper;
                        info_text = 'Improved upper bound';
                        cost = cost-1e-10; % Otherwise we'll fathome!
                    end        
                end
                
                % UPDATE THE LOWER BOUND                
                if isnan(lower)
                    lower = cost;
                end
                if ~isempty(stack)
                    lower =min(cost,min([stack.lower]));
                else
                    lower = min(lower,cost);             
                end
                
                if cost<upper 
                    output = solveupper(p,p_upper,x,options,uppersolver);
                    xu = evaluate_nonlinear(p,output.Primal);
                    upper_residual = resids(p_upper,xu);                   
                    if output.problem == 0 | (all(upper_residual(1:p_upper.K.f)>=-options.bmibnb.eqtol) & all(upper_residual(1+p_upper.K.f:end)>=options.bmibnb.pdtol))
                        this_upper = f+c'*xu+xu'*Q*xu;
                        if (this_upper < (1-1e-5)*upper) & (this_upper < upper - 1e-5)
                            x_min = xu;
                            upper = this_upper;
                            info_text = 'Improved upper bound';
                        end                    
                    end
                else
                    if doplot;ellipplot(diag([200 25]),1,'b',[p.dpos;-p.depth]);drawnow;end
                    info_text = 'Poor lower bound';
                    keep_digging = 0;
                end
            otherwise
            end
        else
            if doplot;ellipplot(diag([200 25]),1,'r',[p.dpos;-p.depth]);drawnow;end
            info_text = 'Infeasible during box-reduction';
            keep_digging = 0;
            cost = inf;
            feasible = 0;          
        end
        solved_nodes = solved_nodes+1;
        
        if ~isempty(stack)
            [stack,lower] = prune(stack,upper,options,solved_nodes); 
        end            
        lower = min(lower,cost); 
        
        % **********************************
        % CONTINUE SPLITTING?
        % **********************************   
        if keep_digging & max(p.ub(p.branch_variables)-p.lb(p.branch_variables))>options.bmibnb.vartol
            spliton = branchvariable(p,options,x);
            bounds  = partition(p,options,spliton,x_min); 
            
            node_1 = savetonode(p,spliton,bounds(1),bounds(2),-1,x,cost);
            node_2 = savetonode(p,spliton,bounds(2),bounds(3),1,x,cost);
            stack = push(stack,node_1);
            stack = push(stack,node_2);        
        end
        
        % Pick and create a suitable node to continue on
        [p,stack] = selectbranch(p,options,stack,x_min,upper);
        
        if isempty(p)
            if ~isinf(upper)
                relgap = 0;
            end
            depth = 0;          
        else      
            relgap = 100*(upper-lower)/(1+abs(upper));
            depth = p.depth;
        end    
        if options.bmibnb.verbose>0
            ws = whos; %Work-space
            Mb = sum([ws(:).bytes])/1024^2; %Megs
            showprogress(sprintf(['%3d ' info_text repmat(' ',1,35-length(info_text)) ' %8.2f%%  %8.4f  %8.4f  %8.4f  %3d  %3d    %5.2fMB    %4.1f%%  '],solved_nodes,relgap,upper,cost,lower,depth,length(stack),Mb,100-vol_reduction*100),options.bmibnb.verbose)
        end                            
        
        absgap = upper-lower;
        % **************************************
        % Continue?
        % **************************************    
        time_ok = cputime-t_start < options.bmibnb.maxtime;
        iter_ok = solved_nodes < options.bmibnb.maxiter;
        any_nodes = ~isempty(p);   
        relgap_too_big = (isinf(lower) | isnan(relgap) | relgap>100*options.bmibnb.relgaptol);
        absgap_too_big = (isinf(lower) | isnan(absgap) | absgap>options.bmibnb.absgaptol);
        taget_not_met = upper>options.bmibnb.target;
        go_on = taget_not_met & time_ok & any_nodes & iter_ok & relgap_too_big & absgap_too_big ;
    end
    if options.bmibnb.verbose>0
        fprintf('******************************************************************************************************************\n')
        if options.bmibnb.verbose;showprogress([num2str2(solved_nodes,3)  ' Finishing.  Cost: ' num2str(upper) ' Gap: ' num2str(relgap) '%'],options.bnb.verbose);end
        fprintf('******************************************************************************************************************\n')
    end
    
    
    
    function stack = push(stackin,p)
    if ~isempty(stackin)
        stack = [p;stackin];
    else
        stack(1)=p;
    end
    
    function [p,stack] = pull(stack,method,x_min,upper);
    if ~isempty(stack)
        switch method
            case 'maxvol'
                for i = 1:length(stack)
                    vol(i) = sum(stack(i).ub(stack(i).branch_variables)-stack(i).lb(stack(i).branch_variables));
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
        p = [];
    end
    
    
    function s = num2str2(x,d,c);
    if nargin==3
        s = num2str(x,c);
    else
        s = num2str(x);
    end
    s = [repmat(' ',1,d-length(s)) s];
    
    
    function res = resids(p,x)
    res= [];
    if p.K.f>0
        res = -abs(p.F_struc(1:p.K.f,:)*[1;x]); 
    end
    if p.K.l>0
        res = [res;p.F_struc(p.K.f+1:p.K.f+p.K.l,:)*[1;x]];
    end
    if (length(p.K.s)>1) | p.K.s>0
        top = startofSOCPCone(p.K)
        for i = 1:length(p.K.s)
            n = p.K.s(i);
            X = p.F_struc(top:top+n^2-1,:)*[1;x];top = top+n^2;
            X = reshape(X,n,n);
            res = [res;min(eig(X))];
        end
    end
    
    res = [res;min([p.ub-x;x-p.lb])];
    
    
    function [stack,lower] = prune(stack,upper,options,solved_nodes)
    % *********************************
    % PRUNE STACK W.R.T NEW UPPER BOUND
    % *********************************
    if ~isempty(stack)
        toolarge = find([stack.lower]>upper*(1-1e-4));
        if ~isempty(toolarge)
            if options.bnb.verbose;showprogress([num2str2(solved_nodes,3) ' Pruned ' num2str(length(toolarge))  '  nodes'],options.bnb.verbose-1);end
            stack(toolarge)=[];
        end
    end
    if ~isempty(stack)
        lower = min([stack.lower]);
    else
        lower = upper;
    end
    
    function pcut = addmcgormick(p)
    
    pcut = p;
    top = 0;
    row = [];
    col = [];
    val = [];
    F_temp = [];
    for i = 1:size(p.nonlins,1)
        z = p.nonlins(i,1);
        x = p.nonlins(i,2);
        y = p.nonlins(i,3);
        x_lb = p.lb(x);
        x_ub = p.ub(x);
        y_lb = p.lb(y);
        y_ub = p.ub(y);
        top = 0;
        row = [];
        col = [];
        val = [];
        
        if x~=y
            row = [1;1;1;1;2;2;2;2;3;3;3;3;4;4;4;4];
            col = [1 ; z+1 ; x+1 ; y+1 ; 1 ; z+1 ; x+1 ; y+1 ; 1 ; z+1 ; x+1 ; y+1 ; 1 ; z+1 ; x+1 ; y+1];
            val = [x_lb*y_lb;1;-y_lb;-x_lb;x_ub*y_ub;1;-y_ub;-x_ub;-x_ub*y_lb;-1;y_lb;x_ub;-x_lb*y_ub;-1;y_ub;x_lb];
            F_temp = [F_temp;sparse(row,col,val,4,size(pcut.F_struc,2))];
        else
            
            nr = 3;
            row = [1;1;1;2;2 ;2; 3; 3; 3];
            col = [1 ;z+1 ;x+1 ;1 ;z+1 ;x+1 ;1 ;z+1 ;x+1];
            val = [-x_ub*x_lb;-1;x_lb+x_ub;x_lb*y_lb;1;-y_lb-x_lb;x_ub*y_ub;1;-y_ub-x_ub];
            
            F_temp = [F_temp;sparse(row,col,val,nr,1+length(p.c))];
        end
        bounds = [x_lb*y_lb x_lb*y_ub x_ub*y_lb x_ub*y_ub];
        if x==y
            pcut.lb(pcut.nonlins(i,1)) = max(pcut.lb(pcut.nonlins(i,1)),max(0,min(bounds)));
        else
            pcut.lb(pcut.nonlins(i,1)) = max(pcut.lb(pcut.nonlins(i,1)),min(bounds));
        end
        pcut.ub(pcut.nonlins(i,1)) = min(pcut.ub(pcut.nonlins(i,1)),max(bounds));
    end
    
    keep = find(~isinf(F_temp(:,1)));
    F_temp = F_temp(keep,:);
    pcut.F_struc = [F_temp;pcut.F_struc];
    pcut.K.l = pcut.K.l+size(F_temp,1);  
    
    function [p,feasible,lower] = lpbmitighten(p,lower,upper,lpsolver)
    
    % Construct problem with only linear terms
    % and add cuts from lower/ upper bounds
    c = p.c;
    p_test = p;
    p_test.K.s = 0;
    p_test.F_struc = p_test.F_struc(1+p_test.K.f:1:p_test.K.l+p_test.K.f,:);
    
    if ~isnan(lower)
        p_test.F_struc = [-(p.lower-abs(p.lower)*0.01) p_test.c';p_test.F_struc];
    end
    if upper < inf
        p_test.F_struc = [upper+abs(upper)*0.01 -p_test.c';p_test.F_struc];
    end
    
    p_test.F_struc = [p_test.lpcuts;p_test.F_struc];
    p_test.K.l = size(p_test.F_struc,1);
    
    % Add cuts for nonlinear terms
    p_test = addmcgormick(p_test);
    
    p_test.F_struc = [p.F_struc(1:1:p.K.f,:);p_test.F_struc];
    
    
    feasible = 1;
    i = 1;
    
    p_test = clean_bounds(p_test);
    
    j = 1;
    n = ceil(max(p.options.bmibnb.lpreduce*length(p_test.linears),1));
    res = zeros(length(p.lb),1);
    for i = 1:size(p.nonlins,1)
        res(p.nonlins(i,2)) = res(p.nonlins(i,2)) + abs( p.x0(p.nonlins(i,1))-p.x0(p.nonlins(i,2)).*p.x0(p.nonlins(i,3)));
        res(p.nonlins(i,3)) = res(p.nonlins(i,3)) + abs( p.x0(p.nonlins(i,1))-p.x0(p.nonlins(i,2)).*p.x0(p.nonlins(i,3)));
    end
    res = res(p.linears);
    [ii,jj] = sort(abs(res));               
    jj = jj(end-n+1:end);    
    
    while feasible & j<=length(jj)
        i = p_test.linears(jj(j));
        if abs(p.ub(i)-p.lb(i)>0.1)
            p_test.c = eyev(length(p_test.c),i); 
            
            output = feval(lpsolver,p_test);
            if output.problem == 0
                if p_test.lb(i) < output.Primal(i)-1e-5
                    p_test.lb(i) = output.Primal(i);  
                    p_test = updateonenonlinearbound(p_test,i);
                end     
                p_test.c = -eyev(length(p_test.c),i);
                output = feval(lpsolver,p_test);
                if output.problem == 0
                    if p_test.ub(i) > output.Primal(i)+1e-5
                        p_test.ub(i) = output.Primal(i);
                        p_test = updateonenonlinearbound(p_test,i);        
                    end
                end
                if output.problem == 1
                    feasible = 0;
                end   
            end
            if output.problem == 1
                feasible = 0;
            end  
            p_test = clean_bounds(p_test);
        end
        j = j + 1;
    end
    p.lb = p_test.lb;
    p.ub = p_test.ub;
    
    
    function p = updateonenonlinearbound(p,changed_var);
    for i = 1:size(p.nonlins,1)
        x = p.nonlins(i,2);
        y = p.nonlins(i,3); 
        if (x==changed_var) | (y==changed_var)
            z = p.nonlins(i,1);        
            x_lb = p.lb(x);
            x_ub = p.ub(x);
            y_lb = p.lb(y);
            y_ub = p.ub(y);
            bounds = [x_lb*y_lb x_lb*y_ub x_ub*y_lb x_ub*y_ub];    
            if x==y        
                p.lb(p.nonlins(i,1)) = max([p.lb(z) 0 min(bounds)]);
                p.ub(p.nonlins(i,1)) = min(p.ub(z),max(bounds));
            else
                p.lb(p.nonlins(i,1)) = max(p.lb(z),min(bounds));
                p.ub(p.nonlins(i,1)) = min(p.ub(z),max(bounds));    
            end 
        end
    end
    
    
    function p = updatenonlinearbounds(p,changed_var,keepbest);
    % if nargin>1
    %     changed_var
    % else
    %     i = 1:size(p.nonlins,1);
    % end
    for i = 1:size(p.nonlins,1)
        z = p.nonlins(i,1);
        x = p.nonlins(i,2);
        y = p.nonlins(i,3);    
        x_lb = p.lb(x);
        x_ub = p.ub(x);
        y_lb = p.lb(y);
        y_ub = p.ub(y);
        bounds = [x_lb*y_lb x_lb*y_ub x_ub*y_lb x_ub*y_ub];    
        if x==y        
            p.lb(p.nonlins(i,1)) = max([p.lb(z) 0 min(bounds)]);
            p.ub(p.nonlins(i,1)) = min(p.ub(z),max(bounds));
        else
            p.lb(p.nonlins(i,1)) = max(p.lb(z),min(bounds));
            p.ub(p.nonlins(i,1)) = min(p.ub(z),max(bounds));    
        end 
    end
    return
    
    if nargin > 1
        
        for i = 1:size(p.nonlins,1)
            z = p.nonlins(i,1);
            x = p.nonlins(i,2);
            y = p.nonlins(i,3);
            if isempty(changed_var) | (x==changed_var) | (y == changed_var) | nargin==3
                bound_x1 = [p.lb(p.nonlins(i,2));p.ub(p.nonlins(i,2))];
                bound_x2 = [p.lb(p.nonlins(i,3));p.ub(p.nonlins(i,3))];
                bounds = [bound_x1(1)*bound_x2(1) bound_x1(1)*bound_x2(2) bound_x1(2)*bound_x2(1) bound_x1(2)*bound_x2(2)];
                if nargin==3
                    if x==y
                        p.lb(p.nonlins(i,1)) = max([p.lb(p.nonlins(i,1)) 0 min(bounds)]);
                        p.ub(p.nonlins(i,1)) = min(p.ub(p.nonlins(i,1)),max(bounds));
                    else
                        p.lb(p.nonlins(i,1)) = max(p.lb(p.nonlins(i,1)),min(bounds));
                        p.ub(p.nonlins(i,1)) = min(p.ub(p.nonlins(i,1)),max(bounds));    
                    end 
                else
                    if x==y
                        p.lb(p.nonlins(i,1)) = max(0,min(bounds));
                        p.ub(p.nonlins(i,1)) = max(bounds);
                    else
                        p.lb(p.nonlins(i,1)) = min(bounds);
                        p.ub(p.nonlins(i,1)) = max(bounds);    
                    end
                end
            end
        end
    else    
        for i = 1:size(p.nonlins,1)
            z = p.nonlins(i,1);
            x = p.nonlins(i,2);
            y = p.nonlins(i,3);
            bound_x1 = [p.lb(p.nonlins(i,2));p.ub(p.nonlins(i,2))];
            bound_x2 = [p.lb(p.nonlins(i,3));p.ub(p.nonlins(i,3))];
            bounds = [bound_x1(1)*bound_x2(1) bound_x1(1)*bound_x2(2) bound_x1(2)*bound_x2(1) bound_x1(2)*bound_x2(2)];
            if x==y
                p.lb(p.nonlins(i,1)) = max( p.lb(p.nonlins(i,1)) ,max(0,min(bounds)));
                p.ub(p.nonlins(i,1)) = min( p.ub(p.nonlins(i,1)) ,max(bounds));
            else
                p.lb(p.nonlins(i,1)) = max(p.lb(p.nonlins(i,1)),min(bounds));
                p.ub(p.nonlins(i,1)) = min(p.ub(p.nonlins(i,1)),max(bounds));    
            end
        end
    end
    
    % *************************************
    % DERIVE LINEAR CUTS FROM SDPs
    % THESE ARE ONLY USED IN BOXREDUCE
    % *************************************
    function p = addsdpcut(p,x)
    if any(p.K.s)
        top = startofSOCPCone(p.K);
        newcuts = 1;
        newF = [];
        for i = 1:length(p.K.s)
            n = p.K.s(i);
            X = p.F_struc(top:top+n^2-1,:)*[1;x];
            X = reshape(X,n,n);
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
            violations = p.lpcuts*[1;x];
            p.lpcuts = p.lpcuts(violations<0.1,:);
            
            if size(p.lpcuts,1)>15*m
                violations = p.lpcuts*[1;x];
                [i,j] = sort(violations);
                p.lpcuts = p.lpcuts(j(1:15*m),:);
                p.lpcuts = p.lpcuts(end-15*m+1:end,:);
            end
        end
    end
    
    
    function spliton = branchvariable(p,options,x)
    
    % Split if box is to narrow
    width = abs(p.ub(p.branch_variables)-p.lb(p.branch_variables));
    if min(width)/max(width) < 0.1
        [i,j] = max(width);
        spliton = p.branch_variables(j);
    else
        %     res = zeros(length(p.lb),1);
        %     for i = 1:size(p.nonlins,1)
        %         res(p.nonlins(i,2)) = res(p.nonlins(i,2)) + abs( x(p.nonlins(i,1))-x(p.nonlins(i,2)).*x(p.nonlins(i,3)));
        %         res(p.nonlins(i,3)) = res(p.nonlins(i,3)) + abs( x(p.nonlins(i,1))-x(p.nonlins(i,2)).*x(p.nonlins(i,3)));
        %     end
        %     
        %      [ii,jj] = sort(abs(res));               
        %      spliton = jj(end);    
        
        res = x(p.nonlins(:,1))-x(p.nonlins(:,2)).*x(p.nonlins(:,3));
        [ii,jj] = sort(abs(res));               
        v1 = p.nonlins(jj(end),2);
        v2 = p.nonlins(jj(end),3);
        
        acc_res1 = sum(abs(res(find((p.nonlins(:,2)==v1) |  p.nonlins(:,3)==v1))));
        acc_res2 = sum(abs(res(find((p.nonlins(:,2)==v2) |  p.nonlins(:,3)==v2))));
        
        if (acc_res1>acc_res2) & ismember(v1,p.branch_variables)
            spliton = v1;
        else
            spliton = v2;
        end         
    end
    
    function bounds = partition(p,options,spliton,x_min)
    
    switch options.bmibnb.branchrule
        case 'omega'
            if ~isempty(x_min)
                bounds = [p.lb(spliton) 0.5*max(p.lb(spliton),min(x_min(spliton),p.ub(spliton)))+0.5*(p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)];                                          
            else
                bounds = [p.lb(spliton) (p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)];                                          
            end
        case 'bisect'
            bounds = [p.lb(spliton) (p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)];                                          
        otherwise
            bounds = [p.lb(spliton) (p.lb(spliton)+p.ub(spliton))/2 p.ub(spliton)]; 
    end
    
    
    
    function [p,feasible,vol_reduction] = boxreduce(p,upper,lower,lpsolver,options);
    
    if options.bmibnb.lpreduce
        
        vol_start    = prod(p.ub(p.branch_variables)-p.lb(p.branch_variables));
        diag_before  =  sum(p.ub(p.branch_variables)-p.lb(p.branch_variables));
        diag_before0 = diag_before;
        
        [pcut,feasible,lower] = lpbmitighten(p,lower,upper,lpsolver);
        diag_after = sum(pcut.ub(p.branch_variables)-pcut.lb(p.branch_variables));      
        iterations = 0;
        while (diag_after/(1e-18+diag_before) < 0.75    ) & feasible & iterations<4
            [pcut,feasible,lower] = lpbmitighten(pcut,lower,upper,lpsolver);
            diag_before = diag_after;
            diag_after = sum(pcut.ub(p.branch_variables)-pcut.lb(p.branch_variables));        
            iterations = iterations + 1;
        end        
        
        % Clean up...
        for i = 1:length(pcut.lb)
            if (pcut.lb(i)>pcut.ub(i)) & (pcut.lb-pcut.ub < 1e-3)
                pcut.lb(i)=pcut.ub(i);
                pcut = updatenonlinearbounds(pcut,i);            
            end
        end
        p.lb = pcut.lb;
        p.ub = pcut.ub;    
        
        % Metric = (V0/V)^(1/n)
        vol_reduction = max(0,min(1,(prod(p.ub(p.branch_variables)-p.lb(p.branch_variables))/(1e-12+vol_start))^(1/length(p.branch_variables))));
        p.lb(p.lb<-1e12) = -inf;
        p.ub(p.ub>1e12) = inf;
    else
        vol_reduction = 1;
        feasible = 1;
    end
    
    function output = solvelower(p,options,lowersolver)
    
    % ********************************************
    % Convex envelope
    % ********************************************
    %p.binary_variables = [];
    p_with_bilinear_cuts = p;
    p_with_bilinear_cuts.F_struc(1:p.K.f,:)=[];
    p_with_bilinear_cuts = addmcgormick(p_with_bilinear_cuts);   
    p_with_bilinear_cuts.F_struc = [p.F_struc(1:p.K.f,:);p_with_bilinear_cuts.F_struc];    
    
    % **************************************
    % SOLVE NODE PROBLEM 
    % **************************************
    if any(p_with_bilinear_cuts.ub<p_with_bilinear_cuts.lb) 
        output.problem=1;
    else
        % We are solving relaxed problem (penbmi might be local solver)
        p_with_bilinear_cuts.monomtable = eye(length(p_with_bilinear_cuts.c));
        
        % fix implied  from mccormick
        % p_with_bilinear_cuts.lb(p.linears) = -inf;
        % p_with_bilinear_cuts.ub(p.linears) = inf;
        % p_with_bilinear_cuts.lb(p.nonlins(:,1)) = -inf;
        % p_with_bilinear_cuts.ub(p.nonlins(:,1)) = inf;
        
        output = feval(lowersolver,p_with_bilinear_cuts); 
        
        relaxed_residual = resids(p_with_bilinear_cuts,output.Primal);
        % Minor clean-up
        output.Primal(output.Primal<p.lb) = p.lb(output.Primal<p.lb);
        output.Primal(output.Primal>p.ub) = p.ub(output.Primal>p.ub);    
    end 
    
    function [p,stack] = selectbranch(p,options,stack,x_min,upper);
    switch options.bmibnb.branchmethod
        case 'maxvol'
            [node,stack] = pull(stack,'maxvol',x_min,upper); 
        case 'best'
            [node,stack] = pull(stack,'best',x_min,upper); 
        otherwise
            [node,stack] = pull(stack,'best',x_min,upper); 
    end
    % Copy node data to p
    if isempty(node)
        p = [];
    else
        p.depth = node.depth;
        p.dpos = node.dpos;
        p.lb = node.lb;
        p.ub = node.ub;
        p.lower = node.lower;
        p.lpcuts = node.lpcuts;
        p.x0 = node.x0;
    end
    
    
    
    function output = solveupper(p,p_original,x,options,uppersolver)
    
    p_upper = p_original;
    
    % Pick an initial point (this can be a bit tricky...)
    % Use relaxed point, shifted towards center of box
    if all(x<=p.ub) & all(x>=p.lb)
        p_upper.x0 = 0.1*x + 0.9*(p.lb+p.ub)/2;
    else
        p_upper.x0 = (p.lb+p.ub)/2;
    end
    % Shift towards interior for variables with unbounded lower or upper
    lbinfbounds = find(isinf(p.lb));
    ubinfbounds = find(isinf(p.ub));
    p_upper.x0(ubinfbounds) = x(ubinfbounds)+0.01;
    p_upper.x0(lbinfbounds) = x(lbinfbounds)-0.01;
    ublbinfbounds = find(isinf(p.lb) & isinf(p.ub));
    p_upper.x0(ublbinfbounds) = x(ublbinfbounds);
    % ...expand the current node just slightly
    p_upper.lb = p.lb;
    p_upper.ub = p.ub;
    p_upper.lb(~isinf(p_original.lb)) = 0.99*p.lb(~isinf(p_original.lb))+p_original.lb(~isinf(p_original.lb))*0.01;
    p_upper.ub(~isinf(p_original.ub)) = 0.99*p.ub(~isinf(p_original.ub))+p_original.ub(~isinf(p_original.ub))*0.01;
    p_upper.lb(isinf(p_original.lb)) = p_upper.lb(isinf(p_original.lb)) - 0.001;
    p_upper.ub(isinf(p_original.ub)) = p_upper.ub(isinf(p_original.ub)) + 0.001;
    p_upper.options.saveduals = 0;
    
    % Solve upper bounding problem
    p_upper.options.warmstart = 1;
    output = feval(uppersolver,p_upper);
    % Project into the box (numerical issue)
    output.Primal(output.Primal<p_upper.lb) = p_upper.lb(output.Primal<p_upper.lb);
    output.Primal(output.Primal>p_upper.ub) = p_upper.ub(output.Primal>p_upper.ub);
    
    
    % This one needs a lot of work
    function p = nonlinear_constraint_propagation(p)
    
    for i = 1:size(p.nonlins,1)
        x = p.nonlins(i,2);
        y = p.nonlins(i,3);
        z = p.nonlins(i,1);
        
        if y==x & p.ub(z)>0
            p.ub(x) = min(p.ub(x),sqrt(p.ub(z)));
            p.lb(x) = max(p.lb(x),-sqrt(p.ub(z)));
        end
        
        if p.lb(y)>0 & p.ub(z)>0 & p.ub(x)>0
            p.ub(x) = min(p.ub(x),p.ub(z)/p.lb(y));
        end
        if p.lb(x)>0 & p.ub(z)>0 & p.ub(y)>0
            p.ub(y) = min(p.ub(y),p.ub(z)/p.lb(x));
        end
        
        if p.ub(y)>0 & p.lb(y)>0 & p.lb(z)>0 
            p.lb(x) = max(p.lb(x),p.lb(z)/p.ub(y));
        end
        if p.ub(x)>0 & p.lb(x)>0 & p.lb(z)>0 
            p.lb(y) = max(p.lb(y),p.lb(z)/p.ub(x));
        end    
    end
    
    
    function vars = decide_branch_variables(p)
    
    if p.options.bmibnb.lowrank==0
        nonlinear = find(~(sum(p.monomtable~=0,2)==1 & sum(p.monomtable,2)==1));
        vars      =  find(sum(abs(full(p.monomtable(nonlinear,:))),1));
    else
        pool1 = p.nonlins(1,2);
        pool2 = p.nonlins(1,3);
        
        for i = 2:size(p.nonlins,1)
            v1 = p.nonlins(i,2);
            v2 = p.nonlins(i,3);
            if v1==v2
                % We are fucked
                pool1 = [pool1 v1];
                pool2 = [pool2 v2];
            else
                if ismember(v1,pool1)
                    pool2 = [pool2 v2];
                elseif ismember(v1,pool2)
                    pool1 = [pool1 v2];
                elseif ismember(v2,pool1)
                    pool2 = [pool2 v1];
                elseif ismember(v2,pool2)
                    pool1 = [pool1 v1];
                else
                    % No member yet
                    pool1 = [pool1 v1];
                    pool2 = [pool2 v2];
                end
            end
        end
        pool1 = unique(pool1);
        pool2 = unique(pool2);
        if isempty(intersect(pool1,pool2))
            if length(pool1)<=length(pool2)
                vars = pool1;
            else
                vars = pool2;
            end
        else
            nonlinear          = find(~(sum(p.monomtable~=0,2)==1 & sum(p.monomtable,2)==1));
            vars =  find(sum(abs(full(p.monomtable(nonlinear,:))),1));  
        end
    end
    
    
    function x = evaluate_nonlinear(p,x);
    x(p.nonlins(:,1)) = x(p.nonlins(:,2)).*x(p.nonlins(:,3));
    
    function p = clean_bounds(p)
    % Fix to improve numerica with integer bounds
    %close = find(1e-6>abs(p.ub - round(p.ub)));
    %p.ub(close) = round(p.ub(close));
    close = 1e-6>abs(p.ub - round(p.ub));
    p.ub(close) = round(p.ub(close));
    
    close = 1e-6>abs(p.lb - round(p.lb));
    p.lb(close) = round(p.lb(close));
    
    p.ub(p.binary_variables) = floor(p.ub(p.binary_variables) + 1e-2);
    %p.lb(p.binary_variables) =  ceil(p.lb(p.binary_variables) - 1e-2);
    %p = updatenonlinearbounds(p);
    
    % Nothing coded to do non-linear propagation
    %p = nonlinear_constraint_propagation(p);
    p.lb(p.lb<-1e12) = -inf;
    p.ub(p.ub>1e12) = inf;
    
    
    
    function node = savetonode(p,spliton,bounds1,bounds2,direction,x,cost);
    
    node.lb = p.lb;
    node.ub = p.ub;
    node.lb(spliton) = bounds1;
    node.ub(spliton) = bounds2;
    if direction == -1
        node.dpos = p.dpos-1/(2^sqrt(p.depth));
    else
        node.dpos = p.dpos+1/(2^sqrt(p.depth));
    end
    node.depth = p.depth+1;
    node.x0 = x;
    node.lpcuts = p.lpcuts;
    node.lower = cost;
    
