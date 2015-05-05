function optimal = bisection(Constraints,Objective,options, tolerance,lower,upper)

% Create an optimizer object which solves for a particular value
% of the sought variable
P = optimizer(Constraints, Objective,options,Objective,recover(depends(Constraints)));

% Initialize the lower bound
if nargin < 5 || isempty(lower)
    lower = 0;
end

% Make sure we actually can solve the lower problem
[sol, flag] = P{lower};
if flag ~=0
    error('Cannot initialize the bisection at lower bound')
end
working_sol = sol;
optimal = lower;

% Now find an upper bound by siply increasing a bound until infeasible
if nargin < 6 || isempty(uper)
    upper = lower + 1;
    [sol, flag] = P{upper};
    while flag ~= 1
        upper = upper*2;
        [sol, flag] = P{upper};
        if upper > 1e6
            error('I suspect there is no upper bound. Testing 1e6 and still feasible...')
        end
    end
end

% Perform bisection
while upper - lower > tolerance
    test = (upper + lower)/2;
    [sol, flag] = P{test};
    if flag == 0
       working_sol = sol;
       optimal = test;
       lower = test;
    else
        upper = test;
    end
end

% Assign computed solution
assign(recover(depends(Constraints)),working_sol);
assign(Objective, optimal);