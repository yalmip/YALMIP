function diagnostic = bisection_core(Constraints,Objective,options)
%BISECTION Core engine

% Note, outer layer added to enable arbitrary sign on ojective.
% Code here initially assumed maximization of scalar, so now we tweak the
% code a bit with swtchedsign to run any case. Fugly.

diagnostic.yalmiptime = 0;
diagnostic.solvertime = 0;
diagnostic.info = '';
diagnostic.problem = 0;

the_sign = 1;
if options.bisection.switchedsign
    the_sign = -1;
end

if length(getvariables(Objective)) > 1
    diagnostic.problem = -4;
    return
end

if  ~(isequal(getbase(Objective),[0 1]))
    diagnostic.problem = -4;
    return
end

if nargin < 3 || isempty(options)
    options = sdpsettings;
end

% Initialize the lower bound
% Typically a good guess
lower = 0;

% Silly upper bound
bestUpper = inf;

% Create an optimizer object which solves feasibility problem for a
% particular value of the sought variable 
x = recover(setdiff(depends(Constraints),depends(Objective)));
if isempty(options) || isequal(options.solver,'')
    % We must figure out suitable solver
    [~,~,~,model] = export(replace(Constraints,Objective,pi));
    options.solver = model.solver.tag;
end
P = optimizer(Constraints,[],options,Objective,x);

% Safe-guard aganst using LMILAB
if ~isempty(strfind(struct(P).model.options.solver,'lmilab'))
    disp('Selected solver LMILAB lacks features making it unsuitable for BISECTION')
    disp('See why here https://yalmip.github.io/solver/lmilab/');
    disp('Select another SDP solver or/and install a better SDP solver');       
    error('Cannot proceed due to poor SDP solver');
    return
end

if options.verbose;
    disp(['Selected solver: ' options.solver]);
    fprintf(['Generating initial bound: ' num2str(lower*the_sign)]);
end

% Make sure we actually can solve the lower problem
solvertime = tic;
[sol, flag] = P{lower};
working_sol = [];
diagnostic.solvertime = diagnostic.solvertime + toc(solvertime);
if flag == 1    
    % This was infeasible, hernce we can use it as an upper bound    
    i = 1;
    while flag
        bestUpper = lower;
        lower = lower - 2^(-4+i);i = i+1;
        if options.verbose;        
            fprintf([' (fail), ' num2str(lower*the_sign)]);
        end
        try
            solvertime = tic;
            [sol, flag] = P{lower};
            diagnostic.solvertime = diagnostic.solvertime + toc(solvertime);
            if lower < -1e6
                  if options.verbose;        
                     fprintf([' (fail). Giving up, never feasible\n']);
                  end
                diagnostic.problem = 21;
                diagnostic.info = yalmiperror(diagnostic.problem,'BISECTION');
                return
            end
        catch            
            diagnostic.problem = 1;
            diagnostic.info = yalmiperror(diagnostic.problem,'BISECTION');
            return
        end
    end
    if options.verbose;        
        fprintf([' (ok).']);
    end
elseif flag == 0    
   working_sol = sol;
else
	diagnostic.problem = flag;
	diagnostic.info = yalmiperror(diagnostic.problem,'BISECTION');
    return
end

v = sol;
optimal = lower;
upper = bestUpper;

if isinf(upper)
    upper = lower+1;
    if options.verbose;        
       fprintf([' (ok), ' num2str(upper*the_sign)]);
    end   
    [sol, flag] = P{upper};
    i = 1;
    while ~flag
        % This was feasible, so we can use it as new lower bound        
        lower = upper;working_sol = sol;optimal = upper;
        % Increase upper bound       
        upper = upper + 2^i;i = i+1;
        try                         
            solvertime = tic;
            if options.verbose;        
                fprintf([' (ok), ' num2str(upper*the_sign)]);
            end   
            [sol, flag] = P{upper};
            diagnostic.solvertime = diagnostic.solvertime + toc(solvertime);
        catch
            upper
            diagnostic.problem = -1;
            diagnostic.info = yalmiperror(diagnostic.problem,'BISECTION');
            return
        end
        if upper > 1e6
            if options.verbose;        
                     fprintf([' (oj). Giving up, always feasible!\n']);
            end
            diagnostic.problem = 2;
            diagnostic.info = yalmiperror(diagnostic.problem,'BISECTION');            
            return
        end
    end
    if options.verbose;        
        fprintf([' (fail).']);
    end   
end

if options.verbose
    fprintf(['\n']);
end

% Perform bisection
iter = 1;
if options.verbose
    disp('Iteration  Lower bound    Test           Upper bound    Gap          Solver status at test');
end
while upper - lower > options.bisection.absgaptol
    test = (upper + lower)/2;
    solvertime = tic;
    [sol, flag] = P{test};
    diagnostic.solvertime = diagnostic.solvertime + toc(solvertime);
    if options.verbose;
        if options.bisection.switchedsign
            L = -upper;
            T = -test;
            U = -lower;
        else
            L = lower;
            T = test;
            U = upper;
        end
        if flag           
           if flag~=1 && ~any(isnan(sol))
                assign(x,sol);assign(Objective,test);res = check(Constraints);
                if min(res) >= -1e-6
                  fprintf(' %4.0f :   %12.5E   %12.5E   %12.5E   %12.5E  %s\n',iter,L,T,U,U-L,[yalmiperror(flag) '(looks ok)'] );   
                  flag = 0;
                else
                  fprintf(' %4.0f :   %12.5E   %12.5E   %12.5E   %12.5E  %s\n',iter,L,T,U,U-L,[yalmiperror(flag) '(assumed infeasible)']);                     
                end
           else               
            fprintf(' %4.0f :   %12.5E   %12.5E   %12.5E   %12.5E  %s\n',iter,L,T, U,U-L,yalmiperror(flag));
           end
        else
            fprintf(' %4.0f :   %12.5E   %12.5E   %12.5E   %12.5E  %s\n',iter,L,T, U,U-L,yalmiperror(flag));
        end
    end
    if flag == 0
       working_sol = sol;
       optimal = test;
       lower = test;
    else
        upper = test;
    end
    iter = iter + 1;
end
if options.bisection.switchedsign
    optimal = -optimal;
end
if isempty(working_sol)
    diagnostic.problem = 1;
else
    % Assign computed solution
    assign(x,working_sol);
    assign(Objective,optimal);
end
if options.verbose  
    if diagnostic.problem==0
        disp(['Bisection terminated successfully with objective ' num2str(optimal)]);
    end
end