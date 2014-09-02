function F = set(F,tag)
% SET Obsolete. Use simple overloading with + instead

if nargin == 2
    F = flatten(F);
    for i = 1:length(F.clauses)
        F.clauses{i}.handle = tag;
    end
end