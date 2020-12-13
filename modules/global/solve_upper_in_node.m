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

[upper_residual,~,feasibility] = constraint_residuals(p_upper,z);
feasible = ~isempty(z) & ~any(isnan(z)) & all(upper_residual(1:p_upper.K.f)>=-p.options.bmibnb.eqtol) & all(upper_residual(1+p_upper.K.f:end)>=-p.options.bmibnb.pdtol);

if feasible
    this_upper = p_upper.f+p_upper.c'*z+z'*p_upper.Q*z;
    if isequal(p_upper.K.f,0) && isequal(p_upper.K.l,0) && isequal(p_upper.K.q,0) && isequal(p_upper.K.s,0) 
        % This is a simple bound-constrained problem. We can try a
        % cross-over and see if it improves the objective
        xx = output.Primal;
        j = find(abs(xx-p_upper.lb)<=p.options.bmibnb.vartol);
        xx(j) = p_upper.lb(j);        
        j = find(abs(xx-p_upper.ub)<=p.options.bmibnb.vartol);
        xx(j) = p_upper.ub(j);
        alt_z = apply_recursive_evaluation(p_upper,xx);
        alt_upper = p_upper.f + p_upper.c'*alt_z + alt_z'*p_upper.Q*alt_z;
        if alt_upper < this_upper
            this_upper = alt_upper;
            z = alt_z;         
        end
    end
    if this_upper < upper
        x_min = z;        
        if length(info_text) == 0
            if isinf(upper)
                info_text = 'Solution found by upper solver';
            else
                info_text = 'Improved solution found by upper solver';
            end
        else
            if isinf(upper)
                info_text = [info_text ' | ' 'Solution found by upper solver'];
            else
                info_text = [info_text ' | ' 'Improved solution found by upper solver'];
            end
         end
        upper = this_upper;
        numglobals = numglobals + 1;
    end
elseif output.problem == 0 && p_upper.K.s(1)>0 && ~p_upper.solver.uppersolver.constraint.inequalities.semidefinite.linear && cutiterations > 0
    % Piggy-back some old code. Generate elementwise cuts from
    % vioalated SDP constraints (if there are any) and add these to the
    % global upper bound model.
    p_u = p_upper;
    p_u.lpcuts = [];
    p_u.cutState = [];
    p_u=createsdpcut(p_u,z);
    p_upper.F_struc = [p_upper.F_struc(1:p_upper.K.f,:);p_u.lpcuts;p_upper.F_struc(1+p_upper.K.f:end,:)];
    p_upper.K.l  = p_upper.K.l + size(p_u.lpcuts,1);
    if ~isempty(p_u.socpcuts)
        p_upper.F_struc = [p_upper.F_struc(1:p_upper.K.f+p_upper.K.l,:);p_u.socpcuts;p_upper.F_struc(1+p_upper.K.f+p_upper.K.l:end,:)];
        if isequal(p_upper.K.q,0)
            p_upper.K.q  = repmat(3,1,size(p_u.socpcuts,1)/3);
        else
            p_upper.K.q  = [repmat(3,1,size(p_u.socpcuts,1)/3) p_upper.K.q];
        end
    end
    
    if ~isempty(p_u.lpcuts)
        n = length(p_u.cutState);        
        if length(info_text) == 0
            info_text = ['Added ' num2str(n) ' cut' plurals(n) ' on SDP cone'];
        else
            info_text = [info_text ' | ' 'Added ' num2str(n) ' cut' plurals(n) ' on SDP cone'];
        end
    end
    
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

function s = plurals(n)
if n>1
    s = 's';
else
    s = '';
end