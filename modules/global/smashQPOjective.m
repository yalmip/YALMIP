function p = smashQPOjective(p,idx)       
p.f = p.f + p.c(idx)'*p.lb(idx);
p.c(idx)=[];
if nnz(p.Q)>0
    p.c = p.c + 2*p.Q(find(~removethese),idx)*p.lb(idx);
    p.f = p.f + p.lb(idx)'*p.Q(idx,idx)*p.lb(idx);
    p.Q(:,find(removethese))=[];
    p.Q(find(removethese),:)=[];
else
    p.Q = spalloc(length(p.c),length(p.c),0);
end