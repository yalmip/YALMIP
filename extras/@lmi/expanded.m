function X = expanded(F,state)

F = flatten(F);
if nargin == 1
    if length(F.LMIid) == 0
        X = [];
    else
        for i = 1:length(F.clauses)
            X(i,1) = F.clauses{i}.expanded;
        end
    end
else
    X = F;
    for i = 1:length(F.clauses)
        X.clauses{i}.expanded = state;
    end
end
