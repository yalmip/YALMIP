function S = groupchanceconstraints(F)
G = {};
S = {};
F = flatten(F);
for i = 1:length(F.clauses)
    if ~isempty(F.clauses{i}.confidencelevel)
        G = addgroup(G,F.clauses{i}.jointprobabilistic);
    end
end
for i = 1:length(G)
    S{i} = [];
end

for i = 1:length(F.clauses)
    for j = 1:length(G)
        if isequal(F.clauses{i}.jointprobabilistic,G{j});
            S{j} = [S{j} i];
        end
    end
end
for i = 1:length(S)  
    s.type = '()';
    s.subs{1} = S{i};
    S{i} = subsref(F,s);
end

function G = addgroup(G,g)
if length(G)==0
    G = {g};
else
    i = 1;
    while i<=length(G)
        if isequal(G{i},g)
            return
        end
        i = i+1;
    end
    G{end+1} = g;
end
