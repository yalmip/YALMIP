function p = preprocess_bilinear_bounds(p)

if ~isempty(p.integer_variables)
    for i = 1:size(p.bilinears,1)
        if ismember(p.bilinears(i,2),p.integer_variables)
            if ismember(p.bilinears(i,3),p.integer_variables)
                p.integer_variables = [p.integer_variables p.bilinears(i,1)];
            end
        end
    end
    if p.K.f > 0
        for i = 1:p.K.f
            if all(ismember(p.F_struc(i,:),[0 1 -1]))
                involved = find(p.F_struc(i,2:end));
                % One variable is linear combination of  integer variables
                if (nnz(ismember(involved,p.integer_variables)) == length(involved)-1) & length(involved)>1
                    p.integer_variables = [p.integer_variables involved];
                end
            end
        end
        p.integer_variables  = unique(p.integer_variables );
    end
end

if ~isempty(p.binary_variables)
    for i = 1:size(p.bilinears,1)
        if ismember(p.bilinears(i,2),p.binary_variables)
            if ismember(p.bilinears(i,3),p.binary_variables)
                p.binary_variables = [p.binary_variables p.bilinears(i,1)];
            end
        end
    end
    for i = 1:p.K.f
        if all(p.F_struc(i,:) == fix(p.F_struc(i,:)))
            involved = find(p.F_struc(i,2:end));
            % One variable is linear combination of binary variables
            if (nnz(ismember(involved,p.binary_variables)) == length(involved)-1) & length(involved)>1
                p.integer_variables = [p.integer_variables involved];
            end
        end
    end
    p.binary_variables = unique(p.binary_variables );
end
