function prob = precalcgpstruct(prob)
prob.b = full(prob.b);
try
    % In 2013a they changed the logic in the output. Tell MATLAB to use old
    % style
    [aa,bb,cc] = unique(prob.map,'legacy');
catch
    % If we have an old version, MATLAB would barf at the legacy flag
    [aa,bb,cc] = unique(prob.map);
end
    
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
