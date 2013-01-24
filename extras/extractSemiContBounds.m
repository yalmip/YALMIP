function [F_struc,K,lb,ub,semicont_variables] = extractSemiContBounds(F_struc,K,lb,ub,semicont_variables);
if ~isempty(semicont_variables)
    % Bounds must be placed in LB/UB
    [lb,ub,cand_rows_eq,cand_rows_lp] = findulb(F_struc,K,lb,ub);
    F_struc(K.f+cand_rows_lp,:)=[];
    F_struc(cand_rows_eq,:)=[];
    K.l = K.l-length(cand_rows_lp);
    K.f = K.f-length(cand_rows_eq);
    redundant = find(lb<=0 & ub>=0);
    semicont_variables = setdiff(semicont_variables,redundant);
end