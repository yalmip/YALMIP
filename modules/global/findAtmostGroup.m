function [g,b] = findAtmostGroup(p,s)
g = [];b = [];
if ~isempty(p.atmost)
    for i = 1:length(p.atmost.groups)
        if all(ismember(s,p.atmost.groups{i}))
            g = p.atmost.groups{i};
            b = p.atmost.bounds(i);
            return
        end
    end
end