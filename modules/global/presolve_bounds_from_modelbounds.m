function model = presolve_bounds_from_modelbounds(model,remove);
if ~isempty(model.F_struc)
    [L,U,cand_rows_eq,cand_rows_lp] = findulb(model.F_struc,model.K);
    model.lb = max([model.lb L],[],2);
    model.ub = min([model.ub U],[],2);    
    model.equalitypresolved = 1;
    if nargin > 1
        if remove
            if ~isempty(cand_rows_lp)
                model.F_struc(model.K.f + cand_rows_lp,:) = [];
                model.K.l = model.K.l - length(cand_rows_lp);
            end
            if ~isempty(cand_rows_eq)
                model.F_struc(cand_rows_eq,:) = [];
                model.K.f = model.K.f - length(cand_rows_eq);
            end            
        end
    end
end


