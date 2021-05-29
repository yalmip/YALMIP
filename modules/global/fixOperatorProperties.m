function p = fixOperatorProperties(p)

for i = 1:length(p.evalMap)
    % y = f(x)
    y = p.evalVariables(i);
    x = p.evalMap{i}.variableIndex;
    xL = p.lb(x);
    xU = p.ub(x);
            
    properties = p.evalMap{i}.properties;
    % Can convexity be fixed for xL <= x <= xL?
    if any(xL<xU) && (isa(properties.convexity,'function_handle') || ~isempty(properties.inflection) || strcmp(properties.convexity,'none'))
        vexity = 'none';
        % User-supplied code to derive convexity
        if isa(properties.convexity,'function_handle')
            vexity = properties.convexity(xL,xU);
        end
        % No such code, or code does not work as expected? Check inflection
        if isequal(vexity,'none') && ~isempty(properties.inflection)
            vexity = DeriveVexityFromInflection(properties,xL,xU);
        end
        if ~isequal(vexity,'none')
            % We have managed to fix convexity
            properties.convexity = vexity;
            % Since convexity is known, we can append a convex hull method
            if isempty(properties.convexhull) && ~isempty(properties.derivative)
                if isequal(vexity,'convex')
                    %f0 = @(x)real(eval([p.evalMap{i}.fcn '(x)']));
                    f0 = p.evalMap{i}.properties.function;
                    f = @(xL,xU)createConvexHullMethodConvex(xL,xU,f0,properties.derivative);
                elseif isequal(vexity,'concave')
                    %f0 = @(x)real(eval([p.evalMap{i}.fcn '(x)']));
                    f0 = p.evalMap{i}.properties.function;
                    f = @(xL,xU)createConvexHullMethodConcave(xL,xU,f0,properties.derivative);
                else
                    f = [];
                end
                properties.convexhull = f;
            end
        end
    end
    if any(xL<xU)
        if isa(properties.monotonicity,'function_handle')
            mono = properties.monotonicity(xL,xU);
            if ~isequal(mono,'none')
                properties.monotonicity = mono;
            end
        elseif ~isempty(properties.stationary) && (isequal(properties.shape,'bell-shape') || isequal(properties.shape,'v-shape'))
            mono = DeriveMonotonicityFromShape(properties,xL,xU);
            if ~isequal(mono,'none')
                properties.monotonicity = mono;
            end
        end
    end
    % Save possibly updated properties
    p.evalMap{i}.properties = properties;
    
    if isempty(p.evalMap{i}.properties.f_upper)
        if isequal(p.evalMap{i}.properties.convexity,'convex')
            %f = makefunction(p.evalMap{i}.fcn);
            f = p.evalMap{i}.properties.function;
            f_upper = @(z,xL,xU)(f(xL) + (z-xL)*(f(xU)-f(xL))/(xU-xL));            
            df_upper = @(z,xL,xU)((f(xU)-f(xL))/(xU-xL));            
            p.evalMap{i}.properties.f_upper = f_upper;
            p.evalMap{i}.properties.df_upper = df_upper;
        end
    end
    if isempty(p.evalMap{i}.properties.f_lower)
         if isequal(p.evalMap{i}.properties.convexity,'concave')
            %f = makefunction(p.evalMap{i}.fcn);
            f = p.evalMap{i}.properties.function;
            f_lower = @(z,xL,xU)(f(xL) + (z-xL)*(f(xU)-f(xL))/(xU-xL));            
            df_lower = @(z,xL,xU)((f(xU)-f(xL))/(xU-xL));            
            p.evalMap{i}.properties.f_lower = f_lower;
            p.evalMap{i}.properties.df_lower = df_lower;
        end
    end
end