function Fout = extractglobalboundsmeta(F)

Fout = [];
F = flatten(F);
meta = find(is(F,'meta'));
Fnew = {};
for i = meta(:)'
    if isequal(F.clauses{i}.data{1},'implies')
        S = F.clauses{i}.data{end};
        if isa(S,'lmi') || isa(S,'constraint')
            Fout = [Fout, extsubsref(S,'Global bound')];
        end
    end
end
