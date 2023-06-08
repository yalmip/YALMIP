function p = presolve_impliedintegers_from_equalities(p)
% Find things like z = 2*x1+3*x2 where x is integer
% and thus implies that z is integer (we only save information
% we do not promote z to integer at currently, hence branching etc
% will only be done on x)
if p.feasible && any(p.K.f)
    implied_integers = [];
    integers = [p.binary_variables p.integer_variables];
    for i = 1:p.K.f
        if fix(p.F_struc(i,1)) == p.F_struc(i,1)
            s = find(p.F_struc(i,2:end));
            if all(s == fix(s))
                y = ismember(s,integers);
                if nnz(y) == length(s)-1
                    newinteger = s(y == 0);
                    implied_integers = [implied_integers newinteger];
                    integers = [integers newinteger];
                end
            end
        end
    end
    p.implied_integers = unique(implied_integers);
else
    p.implied_integers = [];
end