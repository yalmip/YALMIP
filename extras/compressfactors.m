function [Lc,Mc,Rc] = compressfactors(L,M,R)

% Author Johan Löfberg
% $Id: compressfactors.m,v 1.3 2009-11-03 09:36:24 joloef Exp $

Lc = {};
Rc = {};
Mc = {};
taken = zeros(1,length(M));
for i = 1:length(M)
    Lsameas = [];
    Rsameas = [];
    if ~taken(i)
        for j = i+1:length(M)
            if isequal(M{i},M{j})
                if isequal(L{i},L{j})
                    Lsameas = [Lsameas j];
                    taken(j) = 1;
                elseif isequal(R{i},R{j})
                    Rsameas = [Rsameas j];
                    taken(j) = 1;
                end
            end
        end
        if isempty(Rsameas) & ~isempty(Lsameas)
            Lc{end+1} = L{i};
            Mc{end+1} = M{i};
            Rc{end+1} = R{i};
            for j = 1:length(Lsameas)
                Rc{end} = Rc{end} + R{Lsameas(j)};
            end
        elseif isempty(Lsameas) & ~isempty(Rsameas)
            Lc{end+1} = L{i};
            Mc{end+1} = M{i};
            Rc{end+1} = R{i};
            for j = 1:length(Rsameas)
                Lc{end} = Lc{end} + L{Rsameas(j)};
            end
        elseif ~taken(i)
            Lc{end+1} = L{i};
            Mc{end+1} = M{i};
            Rc{end+1} = R{i};
        end
    end
end
