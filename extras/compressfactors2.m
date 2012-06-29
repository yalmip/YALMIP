function [constants,general,newsingles,newpairs] =  compressfactors2(constants,general,singles,pairs)

% Author Johan Löfberg
% $Id: compressfactors2.m,v 1.1 2009-11-03 10:30:33 joloef Exp $

newpairs = {};
taken = zeros(1,length(pairs));
for i = 1:length(pairs)
    Lsameas = [];
    Rsameas = [];
    LsameasR = [];
    RsameasL = [];
    
    if ~taken(i)
        for j = i+1:length(pairs)
            %             if ~isequal(pairs{i}.M,pairs{j}.M)
            %                 if isequal(pairs{i}.M,pairs{j}.M',1)
            %                     pairs{j}.M = pairs{i}.M;
            %                     temp = pairs{j}.L;
            %                     pairs{j}.L = pairs{j}.R';
            %                     pairs{j}.R = temp';
            %                 end
            %             end
            if isequal(pairs{i}.M,pairs{j}.M,1)
                if isequal(pairs{i}.L,pairs{j}.L)
                    Lsameas = [Lsameas j];
                    taken(j) = 1;
                elseif isequal(pairs{i}.R,pairs{j}.R)
                    Rsameas = [Rsameas j];
                    taken(j) = 1;
                 elseif isequal(pairs{i}.M,pairs{i}.M',1)
                     if isequal(pairs{i}.L,pairs{j}.R')
                         temp = pairs{j}.R;
                         pairs{j}.R = pairs{j}.L';
                         pairs{j}.L = temp';
                         Lsameas = [Lsameas j];
                         taken(j) = 1;
                     elseif isequal(pairs{i}.R,pairs{j}.L')
                         temp = pairs{j}.R;
                         pairs{j}.R = pairs{j}.L';
                         pairs{j}.L = temp';
                         Rsameas = [Rsameas j];
                         taken(j) = 1;
                     end
                end
            end
        end
        if isempty(Rsameas) & ~isempty(Lsameas)
            newpairs{end+1}.L = pairs{i}.L;
            newpairs{end}.M = pairs{i}.M;
            newpairs{end}.R = pairs{i}.R;
            for j = 1:length(Lsameas)
                newpairs{end}.R = newpairs{end}.R + pairs{Lsameas(j)}.R;
            end
        elseif isempty(Lsameas) & ~isempty(Rsameas)
            newpairs{end+1}.L = pairs{i}.L;
            newpairs{end}.M = pairs{i}.M;
            newpairs{end}.R = pairs{i}.R;
            for j = 1:length(Rsameas)
                newpairs{end}.L = newpairs{end}.L + pairs{Rsameas(j)}.L;
            end           
        elseif ~taken(i)
            newpairs{end+1}.L = pairs{i}.L;
            newpairs{end}.M = pairs{i}.M;
            newpairs{end}.R = pairs{i}.R;
        end
    end
end

newsingles = singles;
