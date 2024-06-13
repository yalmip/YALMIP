function [p,x_min,upper] = initializesolution(p,x_min,upper)

if p.options.warmstart
    x = p.x0;
    z = evaluate_nonlinear(p,x);
    residual = constraint_residuals(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=-p.options.bmibnb.pdtol);
    if relaxed_feasible 
        % upper_ = p.f+p.c'*z+z'*p.Q*z;
        upper_ = computecost(p.f,p.c,p.Q,z,p);
        if upper_ <= upper
            upper = upper_;
            x_min = z;
        end
    end
else
    % Save for later
    startx0 = p.x0;
    
    p.x0 = zeros(length(p.c),1);
    % Avoid silly warnings
    if ~isempty(p.evalMap)
        for i = 1:length(p.evalMap)
            if (isequal(p.evalMap{i}.fcn,'log') || isequal(p.evalMap{i}.fcn,'log2') || isequal(p.evalMap{i}.fcn,'log10'))
                p.x0(p.evalMap{i}.variableIndex) = (p.lb(p.evalMap{i}.variableIndex) +  p.ub(p.evalMap{i}.variableIndex))/2;
            end
        end
    end
    p = correctEXPConeClosureInitial(p);
    x = p.x0;
    z = evaluate_nonlinear(p,x);
    z = propagateAuxilliary(p,z);
    z = sdpextendsolution(p,z);
    residual = constraint_residuals(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=-p.options.bmibnb.pdtol);
    if relaxed_feasible
        for i = 1:length(p.evalMap)
            if ~isempty(p.evalMap{i}.properties.forbidden)
                if z(p.evalMap{i}.variableIndex) > p.evalMap{i}.properties.forbidden(1) && z(p.evalMap{i}.variableIndex) < p.evalMap{i}.properties.forbidden(2)
                    relaxed_feasible = 0;
                end
            end
        end
    end
    if relaxed_feasible
        infs = isinf(z);
        if ~any(infs)
            % upper_ = p.f+p.c'*z+z'*p.Q*z;
            upper_ = computecost(p.f,p.c,p.Q,z,p);
            if upper_ <= upper
                upper = upper_;
                x_min = z;
                x0 = x_min;
            end
        else
            % Allow inf solutions if variables aren't used in objective
            if all(p.c(infs)==0) & nnz(p.Q(infs,:))==0
                ztemp = z;ztemp(infs)=0;
                %upper_ = p.f+p.c'*ztemp+ztemp'*p.Q*ztemp;
                upper_ = computecost(p.f,p.c,p.Q,ztemp,p);
                if upper_ <= upper
                    upper = upper_;
                    x_min = z;
                    x0 = x_min;
                end
            end
        end
    end
    p.x0 = (p.lb + p.ub)/2;
    if p.binarycardinality.up < length(p.binary_variables)
        % Greedy, set the cheapest ones to 1
        c = p.c(p.binary_variables);
        [~,index] = sort(c,'ascend');
        x0(p.binary_variables(index(1:p.binarycardinality.up))) = 1;
        x0(p.binary_variables(index(p.binarycardinality.up+1:end))) = 0;
    end
    p.x0(isinf(p.ub))=p.lb(isinf(p.ub));
    p.x0(isinf(p.lb))=p.ub(isinf(p.lb));
    p.x0(isinf(p.lb) & isinf(p.ub))=0;
    if ~isempty(p.integer_variables)
        p.x0(p.integer_variables) = round(p.x0(p.integer_variables));
    end
    if ~isempty(p.binary_variables)
        p.x0(p.binary_variables) = round(p.x0(p.binary_variables));
    end
    p = correctEXPConeClosureInitial(p);
    x = p.x0;
    x(isinf(x))=eps;
    x(isnan(x))=eps;
    z = evaluate_nonlinear(p,x);
    z = propagateAuxilliary(p,z);  
    z = sdpextendsolution(p,z);
    residual = constraint_residuals(p,z);
    relaxed_feasible = all(residual(1:p.K.f)>=-p.options.bmibnb.eqtol) & all(residual(1+p.K.f:end)>=p.options.bmibnb.pdtol);   
    if relaxed_feasible 
        tempcost = computecost(p.f,p.c,p.Q,z,p);
        if tempcost < upper
            upper = tempcost;
            x_min = z;
            x0 = x_min;
        end
    end    
    p.x0 = startx0;
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