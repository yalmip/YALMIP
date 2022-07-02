function p = detect_special_objectives(p)
p.LinearBinaryPositiveCost = all(p.c>=0) && all(ismember(find(p.c),p.binary_variables)) && nnz(p.Q)==0 && isempty(p.evalMap);
p.LinearBinaryNegativeCost = all(p.c<=0) && all(ismember(find(p.c),p.binary_variables)) && nnz(p.Q)==0 && isempty(p.evalMap);
p.LinearBinaryCost = all(ismember(find(p.c),p.binary_variables)) && nnz(p.Q)==0 && isempty(p.evalMap);
p.UnitContinuousCost = nnz(p.c) == 1 && all(ismember(find(p.c),p.noninteger_variables)) && nnz(p.Q)==0 && isempty(p.evalMap);
