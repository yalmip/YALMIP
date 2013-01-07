function [upper,x_min,info_text,numglobals,timing] = solve_upper_in_node(p,p_upper,x,upper,x_min,uppersolver,info_text,numglobals,timing);

[output,timing] = global_solve_upper(p,p_upper,x,p.options,uppersolver,timing);
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
        info_text = 'Improved solution';
        numglobals = numglobals + 1;
    end
end