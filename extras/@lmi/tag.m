function F = tag(F,text)
% TAG
%
% Sets the tag on a constraint.

F = flatten(F);
if nargin == 1
    F = F.clauses{1}.handle;
else
for i = 1:length(F.clauses)
    F.clauses{i}.handle = text;
end
end
