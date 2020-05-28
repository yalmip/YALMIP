function p = preprocess_bilinear_bounds(p)

if isempty(p.ub)
    p.ub = repmat(inf,length(p.c),1);
end
if isempty(p.lb)
    p.lb = repmat(-inf,length(p.c),1);
end
if ~isempty(p.F_struc)
    [lb,ub,used_rows_eq,used_rows_lp] = findulb(p.F_struc,p.K);
    if ~isempty([used_rows_eq;used_rows_lp])
        lower_defined = find(~isinf(lb));
        if ~isempty(lower_defined)
            p.lb(lower_defined) = max(p.lb(lower_defined),lb(lower_defined));
        end
        upper_defined = find(~isinf(ub));
        if ~isempty(upper_defined)
            p.ub(upper_defined) = min(p.ub(upper_defined),ub(upper_defined));
        end
        % Remove linear bound inequalities
        if ~isempty(used_rows_lp)
            used_rows_lp = used_rows_lp(find(~any(p.F_struc(p.K.f+used_rows_lp,1+p.nonlinears),2)));
            not_used_rows = setdiff(1:p.K.l,used_rows_lp);
             newKCutl = [];
            for i = 1:length(p.KCut.l)
                 newKCutl  = [newKCutl  find(not_used_rows==p.KCut.l(i))];
               % p.KCut.l(i) = find(not_used_rows == p.KCut.l(i));
               % p.originalModel.KCut.l(i) = find(not_used_rows == p.originalModel.KCut.l(i) );
            end
             p.KCut.l = newKCutl;
            if ~isempty(used_rows_lp)
                p.F_struc(p.K.f+used_rows_lp,:)=[];
              %  p.originalModel.F_struc(p.originalModel.K.f+used_rows_lp,:)=[];
                p.K.l = p.K.l - length(used_rows_lp);
               % p.originalModel.K.l = p.originalModel.K.l - length(used_rows_lp);
            end
        end
        % Remove linear bound inequalities
        if ~isempty(used_rows_eq)
            used_rows_eq = used_rows_eq(find(~any(p.F_struc(used_rows_eq,1+p.nonlinears),2)));                       
            not_used_rows = setdiff(1:p.K.f,used_rows_eq);
            newKCutf = [];
            for i = 1:length(p.KCut.f)
                newKCutf  = [newKCutf  find(not_used_rows==p.KCut.f(i))];
               % p.KCut.f(i) = find(not_used_rows==p.KCut.f(i));
              %  p.originalModel.KCut.f(i) = find(not_used_rows==p.originalModel.KCut.f(i));
            end
            p.KCut.f = newKCutf;
            if ~isempty(used_rows_eq)
                p.F_struc(used_rows_eq,:)=[];
              %  p.originalModel.F_struc(used_rows_eq,:)=[];
                p.K.f = p.K.f - length(used_rows_eq);
              %  p.originalModel.K.f = p.originalModel.K.f - length(used_rows_eq);
            end
        end
    end
end
p.lb(p.binary_variables) = max(0,p.lb(p.binary_variables));
p.ub(p.binary_variables) = min(1,p.ub(p.binary_variables));
p.lb(p.integer_variables) = ceil(p.lb(p.integer_variables));
p.ub(p.integer_variables) = floor(p.ub(p.integer_variables));
p = clean_bounds(p);
if ~isempty(p.bilinears)
    p = propagate_bounds_from_monomials(p);
end