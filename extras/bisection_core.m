function diagnostic = bisection_core(Constraints,Objective,options,tolerance,lower,upper,switchedsign)
%BISECTION Core engine

% Note, outer layer added to enable arbitrary sign on ojective.
% Code here initially assumed maximization of scalar, so now we tweak the
% code a bit with swtchedsign to run any case. Fugly.


if length(getvariables(Objective)) > 1
    error('The objective should be a simple variable');
end

if  ~(isequal(getbase(Objective),[0 1]))
     error('The objective should be a simple  variable');
end

if nargin < 3 || isempty(options)
    options = sdpsettings;
end
% if nargin < 3 || isempty(options)
%     error('An options structure with a specified solver has to be supplied.');
% elseif strcmp(options.solver,'')
%     error('An options structure with a specified solver has to be supplied.');
% end

diagnostic.problem = 0;

if nargin < 4 || isempty(tolerance)
    tolerance = 1e-4;
end

% Initialize the lower bound
if nargin < 5 || isempty(lower)
    lower = 0;
    userGaveLower = 0;
else
    userGaveLower = 1;
end

% Some temp. coding around the fact that we previously assumed maxmization,
% but know use the code for minimization default.
if nargin==7
    switchedSign = 1;   
else
    switchedSign = 1;    
end

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

% Make sure we actually can solve the lower problem
[sol, flag] = P{lower};
if flag ~=0 && userGaveLower
    error('Cannot initialize the bisection at supplied lower bound.')
elseif flag ~=0
    % This was infeasible, hernce we can use it as an upper bound    
    i = 1;
    while flag
        bestUpper = lower;
        lower = lower - 2^i;i = i+1;
        try
            [sol, flag] = P{lower};
            if lower < -1e6
                diagnostic.problem = 1;
                diagnostic.info = yalmiperror(diagnostic.problem,'BISECTION');
            end
        catch            
            diagnostic.problem = -1;
            diagnostic.info = yalmiperror(diagnostic.problem,'BISECTION');
            return
        end
    end
end
v = sol;
optimal = lower;

if nargin == 6 && ~isempty(upper)
    if upper < bestUpper
        % User has supplied an upper bound. Is it really infeasible
        [sol, flag] = P{lower};
        if ~flag
            % No, the supppplied upper bound isn't valid
            upper = bestUpper;
        end
    else
        upper = bestUpper;
    end
else
    upper = bestUpper;
end

if isinf(upper)
    upper = lower+1;
    [sol, flag] = P{upper};
    i = 1;
    while ~flag
        upper = upper + 2^i;i = i+1;
        try             
            [sol, flag] = P{upper};
        catch
            upper
            diagnostic.problem = -1;
            diagnostic.info = yalmiperror(diagnostic.problem,'BISECTION');
            return
        end
        if upper > 1e6
            diagnostic.problem = 2;
            diagnostic.info = yalmiperror(diagnostic.problem,'BISECTION');            
            return
        end
    end
end

% Perform bisection
iter = 1;
if options.verbose;
    disp(['Selected solver: ' options.solver]);
    disp('Iteration  Lower bound    Test           Upper bound    Gap          Status at test');
end
while upper - lower > tolerance    
    test = (upper + lower)/2;
    [sol, flag] = P{test};
    if options.verbose;
        if switchedsign
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
                  fprintf(' %4.0f :   %12.5E   %12.5E   %12.5E   %12.5E  %s\n',iter,L,T,U,U-L,[yalmiperror(flag) '(looks like failure)']);                     
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

% Assign computed solution
assign(x,working_sol);
assign(Objective,optimal);