function F = variablereplace(F,oldVar,newVar)

for i = 1:length(F.clauses)
    if isa(F.clauses{i},'cell')
        for j = 1:length(F.clauses{i})
            F.clauses{i}{j}.data = variablereplace(F.clauses{i}{j}.data,oldVar,newVar);
        end
    else
        F.clauses{i}.data = variablereplace(F.clauses{i}.data,oldVar,newVar);
    end
end