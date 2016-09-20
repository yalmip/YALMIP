function diagnostic = bisection(varargin)
%BISECTION Solve simple quasi-convex MAXIMIZATION problem by bisection
%
%   DIAGNOSTIC = BISECTION(F,h,options,tolerance)
%
%    min   t
%    subject to
%            F(x,t) >=0
%
%   NOTES
%    It is assumed that the problem is quasi-convex in the scalar simple
%    variable t. 
%
%    Lower and upper bounds are automatically detected.
%    Default tolerance 1e-5.
%
%    The algorithm simply solves a series of feasibility problems where the
%    parameter t is fixed, and hones in on optimal value using bisection
%
%    It is recommended to explicitly set a solver. Otherwise YALMIP will
%    have to try to figure ou a suitable solver for the feasibility
%    problems

solvertime = tic;
Constraints = varargin{1};
Objective = varargin{2};
if isequal(getbase(Objective),[0 -1])
    % User wants to maximize something, so we can reuse old code format
    varargin{2} = -Objective;
    options = varargin{3};
    options.bisection.switchedsign = 0;
    options.solver = options.bisection.solver;
    varargin{3} = options;      
    diagnostic = bisection_core(varargin{:});
elseif isequal(getbase(Objective),[0 1])
    % User wants to minimize something. Rewrite as old max code
    options = varargin{3};
    options.bisection.switchedsign = 1;
    options.solver = options.bisection.solver;
    varargin{3} = options;  
    Constraints = replace(Constraints, Objective, -Objective);
    varargin{1} = Constraints;
    diagnostic = bisection_core(varargin{:});   
end
diagnostic.yalmiptime = toc(solvertime)-diagnostic.solvertime;
diagnostic.info = yalmiperror(diagnostic.problem,'bisection');           