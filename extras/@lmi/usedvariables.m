function used = lmi(F)

used = [];
for i = 1:length(F.clauses)
    Fi = F.clauses{i};
    used = unique([used;getvariables(Fi.data)']);
end
