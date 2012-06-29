function [lb,ub,cand_rows] = find_variable_bounds(A,b,Aeq,beq)

n = size(A,2);
lb = -inf*ones(n,1);
ub = inf*ones(n,1);
cand_rows = [];
if size(A,1)>0
    cand_rows = find(sum(A~=0,2)==1);
    for i = 1:length(cand_rows)
        j = find(A(cand_rows(i),:));
        if A(cand_rows(i),j)>0
            ub(j)=min(ub(j),b(cand_rows(i))/A(cand_rows(i),j));
        else
            lb(j)=max(lb(j),b(cand_rows(i))/A(cand_rows(i),j));
        end
    end
end

if nargin>2
    if size(Aeq,2)>0
        candidates = find(sum(Aeq | Aeq,2)==1);
        if length(candidates)>0
            for i = candidates(:)'
                j = find(Aeq(i,:));
                new_bound = beq(i)/Aeq(i,j);
                if any(ub(j) < new_bound) | any(lb(j)> new_bound)
                    infeasible = 1;
                    return
                else
                    ub(j)=new_bound;
                    lb(j)=new_bound;
                end
            end
        end
    end
end