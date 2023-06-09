function d = degree(varargin)

F = varargin{1};
F = flatten(F);
d = subdegree(F.clauses);

function d = subdegree(clauses)
d = -inf;
for i = 1:length(clauses)
    d = max(d,degree(sdpvar(clauses{i}.data)));
end