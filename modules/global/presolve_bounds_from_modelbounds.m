function model = presolve_bounds_from_modelbounds(model,remove)
%% In case we move monomial constraints to bounds
%% we should remember that this was a nonlinear model
model.originallyNonlinearConstraints = 0;
if ~isempty(model.F_struc)      
    if any(any(model.F_struc(:,1+model.nonlinears)))
        model.originallyNonlinearConstraints = 1;
    end
    [L,U,cand_rows_eq,cand_rows_lp] = find_lp_bounds(model.F_struc,model.K);
    model.lb = max([model.lb L],[],2);
    model.ub = min([model.ub U],[],2);    
    model.equalitypresolved = 1;
    if nargin > 1
        if remove
            if ~isempty(cand_rows_lp)
                model.F_struc(model.K.f + cand_rows_lp,:) = [];                
                [~,loc] = ismember(model.KCut.l,setdiff(1:model.K.l,cand_rows_lp));                
                model.KCut.l = loc(find(loc));                                
                model.K.l = model.K.l - length(cand_rows_lp);
            end
            if ~isempty(cand_rows_eq)
                model.F_struc(cand_rows_eq,:) = [];                
                [~,loc] = ismember(model.KCut.f,setdiff(model.KCut.f,cand_rows_eq));
                if ~isempty(loc)
                    model.KCut.f = loc(find(loc));
                end
                model.K.f = model.K.f - length(cand_rows_eq);
            end            
        end
    end
end


