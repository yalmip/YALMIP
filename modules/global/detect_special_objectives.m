function p = detect_special_objectives(p)
notMAXDET = isequal(p.K.m,0) || isempty(p.K.m);
p.LinearBinaryPositiveCost = notMAXDET && all(p.c>=0) && all(ismember(find(p.c),p.binary_variables)) && nnz(p.Q)==0 && isempty(p.evalMap);
p.LinearBinaryNegativeCost = notMAXDET && all(p.c<=0) && all(ismember(find(p.c),p.binary_variables)) && nnz(p.Q)==0 && isempty(p.evalMap);
p.LinearBinaryCost = notMAXDET && all(ismember(find(p.c),p.binary_variables)) && nnz(p.Q)==0 && isempty(p.evalMap);
p.UnitContinuousCost = notMAXDET && nnz(p.c) == 1 && all(ismember(find(p.c),p.noninteger_variables)) && nnz(p.Q)==0 && isempty(p.evalMap);
