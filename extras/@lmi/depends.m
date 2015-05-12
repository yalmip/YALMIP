function LinearVariables = depends(F)

F = flatten(F);

% Get all used variables in this LMI object
used = recursivedepends(F.clauses);

% Now, determine the involved linear variables
[mt,variabletype] = yalmip('monomtable');
if any(variabletype)%1%~isempty(nlv) & any(ismembc(used,nlv(1,:)))
    LinearVariables = find(any(mt(used,:),1));    
    LinearVariables = LinearVariables(:)';
else
    LinearVariables = used;
end

function used = recursivedepends(clauses)

if length(clauses) > 2
    mid = floor(length(clauses)/2);
    used1 = recursivedepends({clauses{1:mid}});
    used2 = recursivedepends({clauses{mid+1:end}});
    used = uniquestripped([used1 used2]);
else
    used = [];
    for i = 1:length(clauses)
        Fivar = getvariables(clauses{i}.data);
        used = uniquestripped([used Fivar(:)']);
    end
end