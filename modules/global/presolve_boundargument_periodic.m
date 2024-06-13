function model = presolve_boundargument_periodic(model)

if ~isempty(model.evalMap)
    % Find arguments only used in periodic functions, and find largest
    % period
    Q = -inf(1,length(model.c));
    for i = 1:length(model.evalMap)
        if ~isempty(model.evalMap{i}.properties.periodic)
            k = model.evalMap{i}.variableIndex;        
            Q(k) = max(Q(k), model.evalMap{i}.properties.periodic);            
        else
            k = model.evalMap{i}.variableIndex;        
            Q(k) = inf;
        end
    end    
    for i = 1:length(model.evalMap)
        if ~isempty(model.evalMap{i}.properties.periodic)
            k = model.evalMap{i}.variableIndex;
            P = Q(k);
            if ~isinf(Q(k))                               
                    if nnz(model.monomtable(:,k)) == 1
                        % It is  ot used in any monomial
                        if isempty(model.F_struc) || nnz(model.F_struc(:,k+1))==0
                            % and not used in any constraint
                            % Hence, since variable is used nowhere except as
                            % argument in this periodic function, we can bound
                            if isinf(model.lb(k)) && isinf(model.ub(k))
                                model.lb(k) = 0;
                                model.ub(k) = P;
                            elseif ~isinf(model.lb(k))
                                model.ub(k) = min(model.ub(k),model.lb(k)+P);
                            elseif ~isinf(model.ub(k))
                                model.lb(k) = max(model.lb(k),model.ub(k)-P);
                        end
                    end
                end
            end
        end
    end
end


