function simplex_model = find_simplex_models(p);

for i = 1:length(p)
    simplex_model(i)= 0;
    if p{i}.K.f == 0
        continue
    elseif any(p{i}.K.q > 0) | any(p{i}.K.s > 0)
        continue
    elseif p{i}.K.f ~= 1
        continue
    else
        aux = p{i};
        b = aux.F_struc(1,1);
        a = aux.F_struc(1,2:end);
        if all(abs(a) == 1)
            b = b/sign(-a(1));
            aux.F_struc(1:p{i}.K.f,:) = [];
            aux.K.f = 0;
            [aux,lower,upper] = find_simple_variable_bounds(aux);
            if all(lower == 0) & aux.K.l == 0 & all((upper == b) | isinf(upper))
                simplex_model(i)=1;
            end
        else
            continue
        end
    end
end