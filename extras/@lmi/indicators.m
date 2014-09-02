function s = indicators(F)

F = flatten(F);
s = F.clauses{1}.extra.indicators;