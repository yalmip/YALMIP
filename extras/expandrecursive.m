function [F_expand,failure,cause] = expandrecursive(variable,F_expand,extendedvariables,monomtable,variabletype,where,level,options,method,extstruct,goal_vexity,allExtStruct,w)

% Some crazy stuff to help out when hand;ling nonocnvex and other
% non-standard problems
global DUDE_ITS_A_GP ALREADY_MODELLED ALREADY_MODELLED_INDEX REMOVE_THESE_IN_THE_END MARKER_VARIABLES OPERATOR_IN_POLYNOM LUbounds CONSTRAINTCUTSTATE

% Assume no failure due to failed model or convexity propagation
cause = '';
failure = 0;

% We only expand this variable if it hasn't already been modelled in
% another constraint, and that model is sufficicent for the current case.
if  ~alreadydone(getvariables(variable),method,goal_vexity)

    % Request model for this variable
    ext_index = find(getvariables(variable) == extendedvariables);

    % [F_graph,properties,arguments,fcn] = model(variable,method,options,extstruct);
    [properties,F_graph,arguments,fcn] = model(variable,method,options,allExtStruct(ext_index),w);

    % If properties has length longer than one, it means that several
    % definitions have been returned. We should now compare the bounds
    % available in the global variable LUbounds, with the domain for these
    % different definitions.
    if length(properties)>1
        all_definitions = properties;
        properties = properties{1};
        % Assume that the first model is the global one
        for i = 2:length(all_definitions)
            if (all_definitions{i}.domain(1)<= LUbounds(getvariables(arguments),1)) & (all_definitions{i}.domain(2)>= LUbounds(getvariables(arguments),2))
                properties = all_definitions{i};
                break
            end      
        end
    else
        if ~isempty(properties)
            properties = properties{1};
        end
    end
    
    % Bit of a hack (CPLEX sucks on some trivial problems, so we have to track
    % these simple models and fix them in the end). Basically, when a
    % constraint such as ismember(x,B) is used, YALMIP introduces a "marker
    % variable" and returns the constraint marker == 1. This constraint and
    % variable should be removed in the end, it is only used to ensure that
    % the operator code is run. CPLEX can fail on this idiot constraints. 
    if ~isempty(properties)
        if isfield(properties,'extra')
            if strcmp(properties.extra,'marker')
                MARKER_VARIABLES = [MARKER_VARIABLES getvariables(variable)];
            end
        end
    end

    % If we are running a graph representation, we have to make sure that
    % we really are satisfying the convexity constraints    
    if isequal(method,'graph')

        if isempty(properties)
            failure = 1;
            exactmodel = 0;
        else            
            failure = ~isequal(properties.convexity,goal_vexity);
            exactmodel = any(strcmpi(properties.model,{'callback','exact','integer'}));
        end

        if failure & ~options.allownonconvex
            failure = 1;
            cause = ['Expected ' goal_vexity ' function in '  where ' at level ' num2str(level)];
            return
        end
       
        if failure & exactmodel
            % This operator is guaranteed to return the associated value
            % (i.e. it is not modelled by graphs) and can thus be used in
            % the model. However, we have to give up convexity propagation
            % since it didn't satisfy the convexity requirements
            failure = 0;
            method = 'exact';
        else
            if failure & options.allowmilp               
                [properties,F_graph,arguments,fcn] = model(variable,'exact',options,allExtStruct(ext_index),w); 
                               
                if ~isempty(properties)
                    properties = properties{1};
                end
   
                if isempty(F_graph)
                    switch fcn
                        case 'sqrt'
                            cause = [cause 'See https://yalmip.github.io/squareroots/'];
                        case 'norm'
                            cause = [cause 'Only 1- and inf-norms can be used in a nonconvex fashion'];
                        otherwise
                            cause = [cause 'The function ' fcn ' does not have an implementation when entering in a ' goal_vexity ' fashion '];
                    end
                    return
                else
                    failure = 0;                   
                    method = properties.model;
                end
            elseif failure
                cause = ['Expected ' 'goal_vexity' ' function in '  where ' at level ' num2str(level)];
            end
        end
    else
        if isempty(F_graph)
            cause = ['Model not available in ' where ' at level ' num2str(level)];
            failure = 1;
            return
        else
            failure = 0;
            method = 'exact';
        end
    end

    % We save variable models for future use
   if failure == 0
   done = save_model_expansion(method,goal_vexity,F_graph,properties);
   if done
       return
   end
   end

   % Now we might have to recurse
   if isa(arguments,'sdpvar')
       if max(size(arguments))>1
           arguments = reshape(arguments,prod(size(arguments)),1);
       else
           try
               if length(getvariables(arguments)) == 1 && isequal(getbase(arguments),[0 1])
               else
                   arguments = recover(getvariables(arguments));
               end
           catch
           end
       end
   end

    [ix,jx,kx] = find(monomtable(getvariables(variable),:));
    if ~isempty(jx) % Bug in 6.1
        if any(kx>1)
            OPERATOR_IN_POLYNOM = [OPERATOR_IN_POLYNOM extendedvariables(jx(find(kx>1)))];
        end
    end

    % Small pre-processing to speed-up large-scale problems (subsref sloooow)
    % with only linear arguments (such as norm(Ax-b) problems)
    if isa(arguments,'sdpvar')
        argvars = getvariables(arguments);
        if max(argvars)>length(variabletype)
            variabletype = yalmip('variabletype');
        end
        do_not_check_nonlinearity = ~any(variabletype(argvars));
        if do_not_check_nonlinearity
            allvariables = argvars;
            fullbasis = getbase(arguments);
            fullbasis = fullbasis(:,2:end);
            fullbasis_transpose = fullbasis';
        end
    else
        do_not_check_nonlinearity = 0;
    end

    j = 1;
    % Ok, here goes the actual recursive code
    while j<=length(arguments) & ~failure

        if length(extendedvariables)==1
            % Optiize for the special case that the whole model only
            % contains 1 single nonlinear operator. If that i the case, we
            % do not need to do recursive stuff. This can save time in some
            % etremely large simple model
            expressionvariables = [];
        else
            if do_not_check_nonlinearity            
                usedvariables = find(fullbasis_transpose(:,j));
                expressionvariables = allvariables(usedvariables);
            else
                expression = arguments(j);
                expressionvariables = unique([depends(expression) getvariables(expression)]);
            end
        end
        index_in_expression = find(ismembcYALMIP(expressionvariables,extendedvariables));

        if ~isempty(index_in_expression)
            for i = index_in_expression
                if ~milpalreadydone(expressionvariables(i))
                    if do_not_check_nonlinearity                       
                        basis = fullbasis(j,usedvariables((i)));
                    else
                        basis = getbasematrix(expression,expressionvariables(i));
                    end
                    go_convex1 = (basis > 0) &  isequal(goal_vexity,'convex') & isequal(properties.monotonicity,'increasing');
                    go_convex2 = (basis <= 0) &  isequal(goal_vexity,'convex') & isequal(properties.monotonicity,'decreasing');
                    go_convex3 = (basis <= 0) &  isequal(goal_vexity,'concave') & isequal(properties.monotonicity,'increasing');
                    go_convex4 = (basis > 0) &  isequal(goal_vexity,'concave') & isequal(properties.monotonicity,'decreasing');
                    go_convex = (go_convex1 | go_convex2 | go_convex3 | go_convex4) & ~strcmpi(method,'integer');
                    
                    go_concave1 = (basis > 0) &  isequal(goal_vexity,'convex') & isequal(properties.monotonicity,'decreasing');
                    go_concave2 = (basis <= 0) &  isequal(goal_vexity,'convex') & isequal(properties.monotonicity,'increasing');
                    go_concave3 = (basis <= 0) &  isequal(goal_vexity,'concave') & isequal(properties.monotonicity,'decreasing');
                    go_concave4 = (basis > 0) &  isequal(goal_vexity,'concave') & isequal(properties.monotonicity,'increasing');
                    go_concave = (go_concave1 | go_concave2 | go_concave3 | go_concave4) & ~strcmpi(method,'integer');
                        
                    if go_convex
                        [F_expand,failure,cause] = expandrecursive(recover(expressionvariables(i)),F_expand,extendedvariables,monomtable,variabletype,where,level+1,options,method,[],'convex',allExtStruct,w);
                    elseif go_concave
                        [F_expand,failure,cause] = expandrecursive(recover(expressionvariables(i)),F_expand,extendedvariables,monomtable,variabletype,where,level+1,options,method,[],'concave',allExtStruct,w);
                    elseif isequal(properties.monotonicity,'exact')
                        [F_expand,failure,cause] = expandrecursive(recover(expressionvariables(i)),F_expand,extendedvariables,monomtable,variabletype,where,level+1,options,method,[],goal_vexity,allExtStruct,w);
                    else
                        if options.allownonconvex%milp
                            [F_expand,failure,cause] = expandrecursive(recover(expressionvariables(i)),F_expand,extendedvariables,monomtable,variabletype,where,level+1,options,'exact',[],'none',allExtStruct,w);
                        else
                            failure = 1;
                            cause = ['Monotonicity required at ' where ' at level ' num2str(level)];
                        end
                    end

                end
            end
        end
        if ~do_not_check_nonlinearity  & ~DUDE_ITS_A_GP & ~options.expandbilinear & ~options.allownonconvex
            if isa(expression,'sdpvar')
                if degree(expression)~=1 &~is(expression,'sigmonial')
                    [Q,c,f,x,info] = quaddecomp(expression);
                    if info
                        failure = 1;
                        cause = ['Polynomial expression  at ' where ' at level ' num2str(level)];
                    else
                        eigv = real(eig(Q));
                        if ~all(diff(sign(eigv))==0)
                            failure = 1;
                            cause = ['Indefinite quadratic in ' where ' at level ' num2str(level)];
                        else
                            fail1 =  isequal(goal_vexity,'convex')  & all(eigv<=0) &  ~isequal(properties.monotonicity,'decreasing');
                            fail2 =  isequal(goal_vexity,'convex')  & all(eigv>0)  &  ~isequal(properties.monotonicity,'increasing');
                            fail3 =  isequal(goal_vexity,'concave') & all(eigv<=0) &  ~isequal(properties.monotonicity,'increasing');
                            fail4 =  isequal(goal_vexity,'concave') & all(eigv>0)  &  ~isequal(properties.monotonicity,'decreasing');

                            if fail1 | fail3
                                failure = 1;
                                cause = ['Concave quadratic encountered in ' where ' at level ' num2str(level)];
                            elseif fail2 | fail4
                                failure = 1;
                                cause = ['Convex quadratic encountered in ' where ' at level ' num2str(level)];
                            end
                        end
                    end
                end
            end
        end
        j = j+1;
    end
    
    % If the variable modelled originates in a constraint which is a cut in
    % the global solver, all defined constraints are also cuts
    if isa(F_graph,'lmi') | isa(F_graph,'constraint') & CONSTRAINTCUTSTATE
        F_graph = setcutflag(lmi(F_graph),CONSTRAINTCUTSTATE);
    end
    
    F_expand = F_expand + F_graph;
end
