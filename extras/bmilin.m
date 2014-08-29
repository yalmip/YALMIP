function diagnostic = bmilin(F,h,options)
%BMILIN Simple BMI solver based on sequential linearizations
%
% diagnostic = bmilin(F,h,options)
%
% EXTREMELY naive implementation of a local BMI solver
% using linearizations and a simple trust-region.
%
% To be precise, it not only "solves" problems with BMIs 
% but also PMIs (polynomial matrix inequalities).
%
% Moreover, the objective may be nonlinear.
%
% The default behaviour of the solver can be
% altered by using the field bmilin in sdpsettings.
%
% The solver should never be called directly by
% the user, but called from solvesdp using the solver tag 'bmilin'
%
% bmilin.trust          - Trust region norm p in ||x-x0||_p < alpha*||x0||_p+beta [2|inf (2)]
% bmilin.alpha          - Trust region parameter [real > 0 (0.5)]
% bmilin.beta           - Trust region parameter [real > 0 (0)]
% bmilin.solver         - Solver for linearized SDP [Any of the solvers above ('')]
% bmilin.maxiterls      - Maximum number of iterations in line search [integer (10)]
% bmilin.maxiter        - Maximum number of iterations [integer (25)]

% Author Johan Löfberg
% $Id: bmilin.m,v 1.3 2005-04-29 08:05:01 joloef Exp $


disp('***********************************************')
disp('')
disp('Warning : This code is still under development')
disp('and should be seen as an example on how you can')
disp('implement your own local BMI solver.')
disp('')
disp('Note also that this solver is slow since it is')
disp('implemented using high-level YALMIP code')
disp('')
disp('***********************************************')

% Recover all used variables
x = recover(depends(F));

% Set up for local solver
verbose = options.verbose;
options.verbose = max(options.verbose-1,0);
options.solver = options.bmilin.solver;

% Initial values hopefully supplied
if options.usex0
    % Initialize to 0 if not initialized
    not_initial = isnan(double(x));
    if any(not_initial)
        assign(x(find(not_initial)),repmat(0,length(find(not_initial)),1));
    end
else
    % No initials, set to zero
    assign(x,repmat(0,length(x),1));
    F_linear = F(find(is(F,'linear')));
    % Find some non-zero by solving for the linear constraints
    if length(F_linear) > 0
        sol = solvesdp(F_linear,linearize(h)+x'*x,options);
        if sol.problem~=0
            diagnostic.solvertime = 0;
            diagnostic.info = yalmiperror(0,'BMILIN');
            diagnostic.problem = sol.problem;
            return
        end
    end
end

% Outer linearization loop
goon = 1;
outer_iter = 0;
alpha = 1;
problem = 0;
while goon
    
    % Save old iterates and old objective function
    x0 = double(x);
    h0 = double(h);
    
    % Linearize
    Flin = linearize(F);
    
    % add a trust region
    switch options.bmilin.trust
    case 1
    case 2
        Flin = Flin + (cone(x-x0,2*alpha*options.bmilin.alpha*norm(x0,2)+options.bmilin.beta));
    case inf
        Flin = Flin + (x-x0 <= options.bmilin.alpha*norm(x0,'inf')+options.bmilin.beta);
        Flin = Flin + (x-x0 >=-options.bmilin.alpha*norm(x0,'inf')+options.bmilin.beta);
    otherwise
    end
    
    % Solve linearized problem
    sol = solvesdp(Flin,linearize(h),options);
    switch sol.problem
    case 0
        % Optimal decision variables for linearized problem
        xnew = double(x);
        
        % Original problem is not guaranteed to be feasible
        % Simple line-search for feasible (and improving) solution
        alpha = 1;
        inner_iter = 0;
        p = checklmi(F);
        while ((min(p)<-1e-5) | (double(h)>h0*1.0001)) & (inner_iter < 15)
            alpha = alpha*0.8;
            assign(x,x0+alpha*(xnew-x0));
            p = checkset(F);
            inner_iter = inner_iter + 1;
        end
        if inner_iter < 300
            outer_iter = outer_iter + 1;
            if verbose > 0
                disp(sprintf('#%d cost : %6.3f  step : %2.3f',outer_iter,double(h),alpha));
            end        
            goon = (outer_iter < options.bmilin.maxiter);
        else
            problem = 4;
            goon = 0;
        end
        
        
    otherwise
        problem = 4;
        goon = 0;
    end
    
end

diagnostic.solvertime = 0;
diagnostic.info = yalmiperror(problem,'BMILIN');
diagnostic.problem = problem;
