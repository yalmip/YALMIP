function model = compressLifted(model)

if isempty(model.aux_variables) || isempty(model.Aeq)
    return
end
    
[used,liftedIndex] = ismember(setdiff(model.aux_variables,model.evalVariables),model.linearindicies);
if isempty(liftedIndex)
    return
end

linearIndex = setdiff(1:length(model.linearindicies),liftedIndex);
lift.d = [];
lift.T = [];
lift.S = [];
if ~isempty(model.Aeq)
    definingLift = find(any(model.Aeq(:,liftedIndex),2));
    lift.d = [lift.d;-model.beq(definingLift)];
    lift.T = [lift.T;model.Aeq(definingLift,linearIndex)];
    lift.S = [lift.S;-model.Aeq(definingLift,liftedIndex)];
    keep = any(lift.S,1);
    if ~all(keep)
        % FIXME: Sort out this case
        return
    end   
    [i,j,k] = find(lift.S);
    lift.T = lift.T(i,:);
    lift.d = lift.d(i);   
    model.Aeq(definingLift,:) = [];
    model.beq(definingLift,:) = [];
    model.Aeq(:,liftedIndex)=[];
    if isempty(model.Aeq)
        model.Aeq = [];
        model.beq = [];
    end
end
if size(lift.S,1) == size(lift.S,2) & length(lift.S) == length(liftedIndex)
    
    % Equalities defining the lifted variables d + Tx = y   
    if ~isempty(model.A)
        Ax = model.A(:,linearIndex);
        Ay = model.A(:,liftedIndex);
        Ax = Ax + Ay*lift.T;
        b = model.b - Ay*lift.d;
        model.A = Ax;
        model.b = b;
    end
    
    lb = model.lb(liftedIndex);model.lb(liftedIndex)=[];
    ub = model.ub(liftedIndex);model.ub(liftedIndex)=[];
    uselb = find(~isinf(lb));
    if ~isempty(uselb)    
        model.A = [model.A;-lift.T(uselb,:)];
        model.b = [model.b;lift.d(uselb) - lb(uselb)];
    end
    useub = find(~isinf(ub));
    if ~isempty(useub)    
        model.A = [model.A;-lift.T(useub,:)];
        model.b = [model.b;lift.d(useub) - lb(useub)];
    end
    model.x0 = model.x0(linearIndex);
    lift.linearIndex = linearIndex;
    lift.liftedIndex = liftedIndex;
    model.lift = lift;   
else
    model.lift = [];
end