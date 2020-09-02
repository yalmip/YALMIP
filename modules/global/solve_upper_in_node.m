function [upper,x_min,info_text,numglobals,timing,p_upper] = solve_upper_in_node(p,p_upper,x,upper,x_min,uppersolver,info_text,numglobals,timing,cutiterations);

if nargin < 10
    cutiterations = 1;
end

% Call nonlinear solver (if there are SDP cones and non-SDP solver, remove
% them and use cut generation later below
[output,timing] = global_solve_upper(p,pruneSDPCone(p_upper),x,p.options,uppersolver,timing);
output.Primal(p_upper.integer_variables) = round(output.Primal(p_upper.integer_variables));
output.Primal(p_upper.binary_variables) = round(output.Primal(p_upper.binary_variables));

z = apply_recursive_evaluation(p_upper,output.Primal);

upper_residual = constraint_residuals(p_upper,z);
feasible = ~isempty(z) & ~any(isnan(z)) & all(upper_residual(1:p_upper.K.f)>=-p.options.bmibnb.eqtol) & all(upper_residual(1+p_upper.K.f:end)>=p.options.bmibnb.pdtol);

if feasible
    this_upper = p_upper.f+p_upper.c'*z+z'*p_upper.Q*z;
    if (this_upper < (1-1e-5)*upper) & (this_upper < upper - 1e-5)
        x_min = z;
        upper = this_upper;
        info_text = 'Improved solution found by upper solver';
        numglobals = numglobals + 1;
    end
elseif p_upper.K.s(1)>0 && ~p_upper.solver.uppersolver.constraint.inequalities.semidefinite.linear
    % Piggy-back some old code. Generate elementwise cuts from
    % vioalated SDP constraints (if there are any) and add these to the
    % global upper bound model.
    p_u = p_upper;
    p_u.lpcuts = [];
    p_u.cutState = [];
    p_u=createsdpcut(p_u,z);
    p_upper.F_struc = [p_upper.F_struc(1:p_upper.K.f,:);p_u.lpcuts;p_upper.F_struc(1+p_upper.K.f:end,:)];
    p_upper.K.l  = p_upper.K.l + size(p_u.lpcuts,1);
    
    % Call recursively if we added cuts, and user wants to run more
    % than one iteration
    if ~isempty(p_u.cutState) > 0 && cutiterations > 1
        [upper,x_min,info_text,numglobals,timing,p_upper] = solve_upper_in_node(p,p_upper,x,upper,x_min,uppersolver,info_text,numglobals,timing,cutiterations-1);
    end    
end

function p = pruneSDPCone(p)
if ~isempty(p.K.s) || p.K.s(1) > 0
    if ~p.solver.uppersolver.constraint.inequalities.semidefinite.linear
        p.F_struc(end-sum(p.K.s.^2)+1:end,:) = [];
        p.K.s = 0;
    end
end