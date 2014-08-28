function F = colon(F,tag)
% COLON Overloaded

% Allows the syntax (x>0):Tag in order to give names/descriptions to
% constraints

F = flatten(F);
for i = 1:length(F.clauses)
    F.clauses{i}.handle = tag;
end
	