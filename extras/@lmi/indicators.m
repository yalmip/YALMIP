function s = indicators(F)

F = flatten(F);
try
    s = F.clauses{1}.extra.indicators;
catch
    s = [];
end