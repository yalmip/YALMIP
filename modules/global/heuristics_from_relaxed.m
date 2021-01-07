function [upper,x_min,cost,info_text,numglobals] = heuristics_from_relaxed(p_upper,x,upper,x_min,cost,numglobals)

x(p_upper.binary_variables) = round(x(p_upper.binary_variables));
x(p_upper.integer_variables) = round(x(p_upper.integer_variables));

z = apply_recursive_evaluation(p_upper,x(1:length(p_upper.c)));

relaxed_residual = constraint_residuals(p_upper,z);

eq_ok = all(relaxed_residual(1:p_upper.K.f)>=-p_upper.options.bmibnb.eqtol);
iq_ok = all(relaxed_residual(1+p_upper.K.f:end)>=-p_upper.options.bmibnb.pdtol);

relaxed_feasible = eq_ok & iq_ok;
info_text = '';
if relaxed_feasible
    this_upper = p_upper.f+p_upper.c'*z+z'*p_upper.Q*z;
    if (this_upper < (1-1e-5)*upper) & (this_upper < upper - 1e-5)
        x_min = x;
        upper = this_upper;
        info_text = 'Improved solution found by heuristics';
        cost = cost-1e-10; % Otherwise we'll fathome!
        numglobals = numglobals + 1;
    end
end