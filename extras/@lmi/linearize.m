function [F,changed] = linearize(F)
%LINEARIZE Linearizes all constraints

changed = 0;
F = flatten(F);
Counter = length(F.LMIid);
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
