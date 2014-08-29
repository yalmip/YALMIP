function [lb,ub,cand_rows_eq,cand_rows_lp] = findlinearulb(F_struc,K,lb,ub,linearvariables)
%FINDULB Internal function to extract upper and lower variable bounds

n = size(F_struc,2)-1;
if nargin < 3
    lb = -inf*ones(n,1);
    ub = inf*ones(n,1);
end
cand_rows_eq = [];
cand_rows_lp = [];
ub2 = ub;
lb2 = lb;
if (K.f ~=0)
    A = -F_struc(1:K.f,2:end);
    b = F_struc(1:K.f,1);
    n = size(F_struc,2)-1;
    cand_rows_eq = find((sum(A~=0,2))==1 & (sum(A(:,linearvariables)~=0,2)==1));
    for i = 1:length(cand_rows_eq)
        j = find(A(cand_rows_eq(i),:));        
        ub(j)=min(ub(j),b(cand_rows_eq(i))/A(cand_rows_eq(i),j));
        lb(j)=max(lb(j),b(cand_rows_eq(i))/A(cand_rows_eq(i),j));
    end
end

if (K.l ~=0)
    A = -F_struc(K.f+1:K.f+K.l,2:end);
    b = F_struc(K.f+1:K.f+K.l,1);
    n = size(F_struc,2)-1;
    cand_rows_lp = find((sum(A~=0,2))==1 & (sum(A(:,linearvariables)~=0,2)==1));
    if ~isempty(cand_rows_lp)
        [ii,jj,kk] = find(A(cand_rows_lp,:));
        s_pos = find(kk>0);
        s_neg = find(kk<=0);
        if ~isempty(s_pos)
            for s = 1:length(s_pos)
                ub(jj(s_pos(s)),1) = full(min(ub(jj(s_pos(s))),b(cand_rows_lp(ii(s_pos(s))))./kk(s_pos(s))));
            end
        end
        if ~isempty(s_neg)
            for s = 1:length(s_neg)
                lb(jj(s_neg(s)),1) = full(max(lb(jj(s_neg(s))),b(cand_rows_lp(ii(s_neg(s))))./kk(s_neg(s))));
            end
        end
    end
end