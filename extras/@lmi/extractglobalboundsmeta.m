function Fout = extractglobalboundsmeta(F)

Fout = [];
F = flatten(F);
meta = find(is(F,'meta'));
Fnew = {};
for i = meta(:)'
    if isequal(F.clauses{i}.data{1},'implies')
        S = F.clauses{i}.data{end};
        if isa(S,'lmi')  
            try
                Fout = [Fout, extsubsref(S,'Global bound')];
            catch
            end
        elseif isa(S,'constraint')
            g = struct(S).tag;
            if isequal(g,{'Global bound'})
                Fout = [Fout, S];
            end
        end
    end
end
