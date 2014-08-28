function LIST = constraintclass(F,property)
%CONSTRAINTCLASS   Returns a list describing the constraints

% Author Johan Löfberg 
% $Id: constraintclass.m,v 1.2 2005-02-04 10:10:26 johanl Exp $   

F = flatten(F);
if isempty(F.clauses)
    LIST = [];   
else
    LIST = [];
    for i = 1:length(F.clauses)
        LIST = [LIST;F.clauses{i}.type];
    end    
end

