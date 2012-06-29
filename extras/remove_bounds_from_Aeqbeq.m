function [lb,ub,A,b,infeasible] = remove_bounds_from_Aeqbeq(A,b,lbin,ubin);

infeasible = 0;
lb = lbin;
ub = ubin;
if size(A,1)>0
    cand_rows_lp = find(sum(A~=0,2)==1);
    if ~isempty(cand_rows_lp)
        [ii,jj,kk] = find(A(cand_rows_lp,:));
        s_pos = find(kk>0);
        s_neg = find(kk<=0);
        if ~isempty(s_pos)
            for s = 1:length(s_pos)
                lb(jj(s_pos(s)),1) = full(b(cand_rows_lp(ii(s_pos(s))))./kk(s_pos(s)));
                ub(jj(s_pos(s)),1) = full(b(cand_rows_lp(ii(s_pos(s))))./kk(s_pos(s)));
            end
        end
        if ~isempty(s_neg)
            for s = 1:length(s_neg)
                lb(jj(s_neg(s)),1) = full(b(cand_rows_lp(ii(s_neg(s))))./kk(s_neg(s)));
                ub(jj(s_neg(s)),1) = full(b(cand_rows_lp(ii(s_neg(s))))./kk(s_neg(s)));
            end
        end
    end
    A(cand_rows_lp,:) = [];
    b(cand_rows_lp,:) = [];
end
if any(lb > ubin)
    j = find(lb > ubin);
    ub(j) = ubin(j);
    infeasible = 1;
end
if any(ub < lbin)
    j = find(ub < lbin);
    lb(j) = lbin(j);
    infeasible = 1;
end
