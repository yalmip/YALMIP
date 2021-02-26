function [p,x_min,upper] = initializesolution(p);

x_min = zeros(length(p.c),1);
upper = inf;
if p.options.usex0
    x = p.x0;
    z = evaluate_nonlinear(p,x);
    residual = constraint_residuals(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=p.options.bmibnb.pdtol);
    if relaxed_feasible
        upper = p.f+p.c'*z+z'*p.Q*z;
        x_min = z;
    end
else
    x0 = p.x0;
    p.x0 = zeros(length(p.c),1);
    % Avoid silly warnings
    if ~isempty(p.evalMap)
        for i = 1:length(p.evalMap)
            if (isequal(p.evalMap{i}.fcn,'log') | isequal(p.evalMap{i}.fcn,'log2') | isequal(p.evalMap{i}.fcn,'log10'))
                p.x0(p.evalMap{i}.variableIndex) = (p.lb(p.evalMap{i}.variableIndex) +  p.ub(p.evalMap{i}.variableIndex))/2;
            end
        end
    end
    x = p.x0;
    z = evaluate_nonlinear(p,x);
    z = propagateAuxilliary(p,z);
    
    residual = constraint_residuals(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=p.options.bmibnb.pdtol);
    if relaxed_feasible
        infs = isinf(z);
        if isempty(infs)
            upper = p.f+p.c'*z+z'*p.Q*z;
            x_min = z;
            x0 = x_min;
        else
            % Allow inf solutions if variables aren't used in objective
            if all(p.c(infs)==0) & nnz(p.Q(infs,:))==0
                ztemp = z;ztemp(infs)=0;
                upper = p.f+p.c'*ztemp+ztemp'*p.Q*ztemp;
                x_min = z;
                x0 = x_min;
            end
        end
    end
    p.x0 = (p.lb + p.ub)/2;
    p.x0(isinf(p.ub))=p.lb(isinf(p.ub));
    p.x0(isinf(p.lb))=p.ub(isinf(p.lb));
    p.x0(isinf(p.lb) & isinf(p.ub))=0;
    if ~isempty(p.integer_variables)
        p.x0(p.integer_variables) = round(p.x0(p.integer_variables));
    end
    if ~isempty(p.binary_variables)
        p.x0(p.binary_variables) = round(p.x0(p.binary_variables));
    end
    
    x = p.x0;
    x(isinf(x))=eps;
    x(isnan(x))=eps;
    z = evaluate_nonlinear(p,x);
    z = propagateAuxilliary(p,z);
    
    residual = constraint_residuals(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=p.options.bmibnb.pdtol);
    if relaxed_feasible & ( p.f+p.c'*z+z'*p.Q*z < upper)
        upper = p.f+p.c'*z+z'*p.Q*z;
        x_min = z;
        x0 = x_min;
    end
    p.x0 = x0;
end


function z = propagateAuxilliary(p,z)

try
    % New feature. If we introduce new variables xx = f(x) to be used
    % in a nonlinear operator, we can derive its value when x is chosen
    
    if ~isempty(p.aux_variables)
        if p.K.f > 1
            A = p.F_struc(1:p.K.f,2:end);
            b = p.F_struc(1:p.K.f,1);
            for i = 1:length(p.aux_variables)
                j = find(A(:,p.aux_variables(i)));
                if length(j)==1
                    if A(j,p.aux_variables(i))==1
                        z(p.aux_variables(i)) = -b(j)-A(j,:)*z;
                    end
                end
            end
            z = evaluate_nonlinear(p,z);
        end
    end
catch
end
