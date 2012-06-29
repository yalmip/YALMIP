function F = eliminateBinary(F,binaries)

for i = 1:length(F.clauses)
    data = eliminateBinary(F.clauses{i}.data,binaries);
    F.clauses{i}.data = data;
end