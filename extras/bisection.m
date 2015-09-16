function optimal = bisection(Constraints,Objective,options,tolerance,lower,upper)
%BISECTION Solve simple quasi-convex MAXIMIZATION problem by bisection
%
%   OPTIMALH = BISECTION(F,h,options,tolerance,lower,upper)
%
%    max   h
%    subject to
%            F(x,h) >=0
%
%   NOTES
%    It is assumed that the problem is quasi-convex in the scalar simple
%    variable h. 
%
%    Lower and upper bounds are automatically detected if not assigned.
%    Default tolerance 1e-5.
%
%    A suitable solver has to be specified in the solver options.


if length(getvariables(Objective)) > 1
    error('The objective should be a simple variable');
end

if  ~(isequal(getbase(Objective),[0 1]))
     error('The objective should be a simple  variable');
end

if nargin < 3 || isempty(options)
    error('An options structure with a specified solver has to be supplied.');
elseif strcmp(options.solver,'')
    error('An options structure with a specified solver has to be supplied.');
end

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

bestUpper = inf;

% Create an optimizer object which solves feasibility problem for a
% particular value of the sought variable 
x = recover(setdiff(depends(Constraints),depends(Objective)));
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
        [sol, flag] = P{lower};
        if lower < -1e6
            error('I suspect there is no lower bound. I cannot find any feasible solution while decreasing lower bound')
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
        [sol, flag] = P{upper};
        if upper > 1e6
            error('I suspect there is no upper bound. I cannot find any infeasible solution while increasing upper bound')
        end
    end
end

% Perform bisection
iter = 1;
if options.verbose;
    disp('Iteration  Lower bound    Test           Upper bound    Gap          Status at test');
end
while upper - lower > tolerance    
    test = (upper + lower)/2;
    [sol, flag] = P{test};
    if options.verbose;
        if flag
           % fprintf(' %4.0f : %12.3E  %12.3E  %12.3E    %s\n',iter,lower,test, upper,'Infeasible');
           if flag~=1 && ~any(isnan(sol))
                assign(x,sol);assign(Objective,test);res = check(Constraints);
                if min(res) >= -1e-6
                  fprintf(' %4.0f :   %12.5E   %12.5E   %12.5E   %12.5E  %s\n',iter,lower,test, upper,upper-lower,[yalmiperror(flag) '(looks ok)'] );   
                  flag = 0;
                else
                  fprintf(' %4.0f :   %12.5E   %12.5E   %12.5E   %12.5E  %s\n',iter,lower,test, upper,upper-lower,[yalmiperror(flag) '(looks like failure)']);                     
                end
           else               
            fprintf(' %4.0f :   %12.5E   %12.5E   %12.5E   %12.5E  %s\n',iter,lower,test, upper,upper-lower,yalmiperror(flag));
           end
        else
            fprintf(' %4.0f :   %12.5E   %12.5E   %12.5E   %12.5E  %s\n',iter,lower,test, upper,upper-lower,yalmiperror(flag));
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