function S = extractRandomDefinitions(F)
F = flatten(F);
S = {};
for i = 1:length(F)
    S{i}.variables = F.clauses{i}.data;
    S{i}.distribution =  struct(F.clauses{i}.data).extra.distribution;
end
