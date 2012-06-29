function diagnostic = callbmialt(F,h,options)
%
% EXTREMELY naive implementation of a local BMI solver
% using alternating directions.
%
% The default behaviour of the solver can be
% altered be using the field bmialt in sdpsettings.
%

% Recover all used variables
      
x_vars   = depends(F);    
all_vars = getvariables(F);

% Nonlinear variables
nonlin_variables = setdiff(all_vars,x_vars);

x_nonlin = recover(nonlin_variables);

set1_vars = [];
set2_vars = [];
for i = 1:length(x_nonlin);
    xi_vars = depends(x_nonlin(i));
    if length(xi_vars)~=2
        error('Not bilinear')
    end
    if isempty(set1_vars)
        set1_vars = xi_vars(1);
        set2_vars = xi_vars(2);
    else
 %       index = [ismember(xi_vars(1),set1_vars) ismember(xi_vars(2),set1_vars) ismember(xi_vars(1),set2_vars) ismember(xi_vars(2),set2_vars)];
    
        if ismember(xi_vars(1),set1_vars)
            set2_vars = [set2_vars xi_vars(2)];
        elseif ismember(xi_vars(1),set2_vars)
            set1_vars = [set1_vars xi_vars(2)];
        elseif ismember(xi_vars(2),set1_vars)
            set2_vars = [set2_vars xi_vars(1)];
        elseif ismember(xi_vars(2),set2_vars)  
            set1_vars = [set1_vars xi_vars(1)];            
        end       
    end
end
    

set1_vars = unique(set1_vars);
set2_vars = unique(set2_vars);

x = recover(x_vars);    
x1 = recover(set1_vars);    
x2 = recover(set2_vars);    

% Set up for local solver
verbose = options.verbose;
options.verbose = max(options.verbose-1,0);
options.solver = options.bmilin.solver;

% Initial values hopefully supplied
if options.usex0
    % Initialize to 0 if not initialized
    not_initial = isnan(double(x));
    if any(not_initial)
        setsdpvar(x(find(not_initial)),repmat(0,length(find(not_initial)),1));
    end
else
    % No initials, set to zero
    setsdpvar(x,repmat(0,length(x),1));
    F_linear = F(find(is(F,'linear')));
    % Find some non-zero by solving for the linear constraints
    if length(F_linear) > 0
        solvesdp(F_linear,linearize(h)+x'*x,options);
    end
end


% Outer linearization loop
goon = 1;
outer_iter = 0;
while goon
    
    % Save old iterates and old objective function
    x0 = double(x);
    h0 = double(h);
    
    Flin0 = replace(F,x1,double(x1));
    obj1 = replace(h,x1,double(x1));
    
    Flin = set([]);
    for i = 1:length(Flin0)
        if isa(sdpvar(Flin0(i)),'sdpvar')
            Flin = Flin + Flin0(i);
        end
    end
    
    % Solve linearized problem
    solvesdp(Flin,obj1,options);

        if verbose > 0
        disp(sprintf('#%d cost : %6.3f  ',outer_iter,double(h)));
    end

    Flin0 = replace(F,x2,double(x2));
    obj1 = replace(h,x2,double(x2));
    
    Flin = set([]);
    for i = 1:length(Flin0)
        if isa(sdpvar(Flin0(i)),'sdpvar')
            Flin = Flin + Flin0(i);
        end
    end
    
    
    % Solve linearized problem
    solvesdp(Flin,obj1,options);

    outer_iter = outer_iter + 1;
    if verbose > 0
        disp(sprintf('#%d cost : %6.3f  ',outer_iter,double(h)));
    end
    goon = (outer_iter < options.bmilin.maxiter);
end

diagnostic.solvertime = 0;
diagnostic.info = yalmiperror(0,'BMILIN');
diagnostic.problem = 0;
