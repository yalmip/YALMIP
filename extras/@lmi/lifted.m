function X = expanded(F,state)

if nargin == 1
    if length(F.clauses) == 0
        X = [];
    else
        for i = 1:length(F.clauses)
            X(i,1) = F.clauses{i}.lift;
        end
    end
else
    X = F;
    for i = 1:length(F.clauses)
        X.clauses{i}.lift = state;
    end
end
