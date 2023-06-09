function p = smashQPOjective(p,removethese)       
p.f = p.f + p.c(removethese)'*p.lb(removethese);
p.c(removethese)=[];
if nnz(p.Q)>0
    p.c = p.c + 2*p.Q(find(~removethese),removethese)*p.lb(removethese);
    p.f = p.f + p.lb(removethese)'*p.Q(removethese,removethese)*p.lb(removethese);
    p.Q(:,find(removethese))=[];
    p.Q(find(removethese),:)=[];
else
    p.Q = spalloc(length(p.c),length(p.c),0);
end