function prob = precalcgpstruct(prob)
prob.b = full(prob.b);
[aa,bb,cc] = unique(prob.map);
if aa(1)==1
    bb=[0 bb(1:end)'];
end
indsi = [];
indsj = [];
vals = [];
for i = 1:max(prob.map)
    ind = [(bb(i)+1):(bb(i+1))];
    indsj = [indsj ind];
end
indsi = prob.map(prob.map>0);
vals = prob.b(prob.map>0);
prob.B = sparse(indsi,indsj,vals,max(prob.map),size(prob.A,1));

ind = find(prob.map==0);
prob.Afun = prob.A(ind,:);
prob.bfun = prob.b(ind);
