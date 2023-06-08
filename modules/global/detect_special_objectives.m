function p = detect_special_objectives(p)
notMAXDET = isequal(p.K.m,0) || isempty(p.K.m);
notQP = nnz(p.Q)==0;
notNONLINEAR = notQP && notMAXDET && isempty(p.evalMap);
binaryCost = all(ismember(find(p.c),p.binary_variables));
p.LinearBinaryPositiveCost = notNONLINEAR && all(p.c>=0) && binaryCost;
p.LinearBinaryNegativeCost = notNONLINEAR && all(p.c<=0) && binaryCost;
p.LinearBinaryCost = notNONLINEAR && binaryCost;
p.UnitContinuousCost = notNONLINEAR && nnz(p.c) == 1 && all(ismember(find(p.c),p.noninteger_variables)) && (length(p.noninteger_variables)==1);

if nnz(p.Q-fix(p.Q))==0 && (nnz(p.c-fix(p.c))==0) && notMAXDET
    discrete_variables = union(p.binary_variables,p.integer_variables);
    can_use_ceil_lower_ = all(ismember(find(p.c),discrete_variables));
    can_use_ceil_lower_ = can_use_ceil_lower_ && all(ismember(find(any(p.Q,2)),discrete_variables));
    p.IntegerCost = can_use_ceil_lower_;
else    
    p.IntegerCost = 0;
end