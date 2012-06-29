function [F,changed] = linearize(F)
%LINEARIZE Linearizes all constraints

% Author Johan Löfberg
% $Id: linearize.m,v 1.3 2005-02-04 10:10:27 johanl Exp $

changed = 0;
Counter = size(F.clauses,2);
for i = 1:Counter
    switch F.clauses{i}.type
        case {1,2,3}
            Fi = F.clauses{i}.data;
            if ~is(Fi,'linear')
                Flin = linearize(Fi);
                F.clauses{i}.data = Flin;
            end
        otherwise
    end
end
