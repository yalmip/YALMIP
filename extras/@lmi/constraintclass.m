function LIST = constraintclass(F,property)
%CONSTRAINTCLASS   Returns a list describing the constraints

F = flatten(F);
if isempty(F.clauses)
    LIST = [];   
else
    LIST = [];
    for i = 1:length(F.clauses)
        LIST = [LIST;F.clauses{i}.type];
    end    
end

