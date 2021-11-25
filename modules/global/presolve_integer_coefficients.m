function p = presolve_integer_coefficients(p)

% We might use these shortuts later
p.isinteger = zeros(length(p.c),1);
p.isinteger(p.binary_variables) = 1;
p.isinteger(p.integer_variables) = 1;

% Trivially stupid equalities (i.e. Jeroslows example...)
if any(p.K.f)
    for i = 1:p.K.f
        r = find(p.F_struc(i,2:end));
        if any(r) && all(p.isinteger(r))            
            a = p.F_struc(i,1+r);
            if all(fix(a)==a)
                % Integer coefficients row               
                m = abs(gcdfactor(a));
                if m~=1
                    p.F_struc(i,:) = p.F_struc(i,:)/m;
                end
                if p.F_struc(i,1) ~= fix(p.F_struc(i,1))
                    p.feasible = 0;
                    return
                end
            end
        end
    end
end
if any(p.K.l)
    top = startofLPCone(p.K);
    for i = 1:p.K.l
        r = find(p.F_struc(top,2:end));
        if any(r) && all(p.isinteger(r))
            a = p.F_struc(top,1+r);
            if all(fix(a)==a)
                % Integer coefficients row
                m = abs(gcdfactor(a));
                if m~=1
                    p.F_struc(top,:) = p.F_struc(top,:)/m;
                end
                p.F_struc(top,1) = floor(p.F_struc(top,1)+1e-10);
            end
        end
        top = top + 1;
    end
end