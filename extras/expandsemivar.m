function  [F,semicont_variables,binary_variables] = expandsemivar(F,semicont_variables,binary_variables)
if isempty(semicont_variables)
    return
else
    LU = getbounds(F);    
    L = LU(semicont_variables,1);
    U = LU(semicont_variables,2);
    if any(isinf(L))
        error('There are semi-continuous variables without explicit (when non-zero) lower bound');
    elseif any(isinf(U))
        error('There are semi-continuous variables without explicit (when non-zero) upper bound');
    end
    n = length(semicont_variables);
    r = binvar(n,1);binary_variables = [binary_variables,getvariables(r)];
    x = recover(semicont_variables);
    F = [F, r.*L <= x <= U.*r];
end