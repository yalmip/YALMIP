function F = expandmeta(F)

F = flatten(F);
meta = find(is(F,'meta'));
Fnew = {};
for i = meta(:)'   
    Fnew{i} = feval(F.clauses{i}.data{1},'expand',[],F.clauses{i}.data{2:end});    
end
F.clauses = {F.clauses{setdiff(1:size(F.clauses,2),meta)}};
F.LMIid(meta) = [];
F = [F,Fnew{:}];