function F = lowrank(F,x)

identifiers = [];
for i = 1:size(F.clauses,2)
    switch F.clauses{i}.type
        case 1
            identifiers = [identifiers F.LMIid(i)];
        otherwise            
    end
end

if isempty(identifiers)
    F = set([]);
    return
end

if nargin>1
    variables = getvariables(x);
else
    variables = [];
end

lrData.id = identifiers;
lrData.variables = variables;

%F.clauses = {{F.clauses{1}}};
F.clauses{1}.data = lrData;
F.clauses{1}.type  = 14;
F.LMIid = -1;
F.clauses{1}.strict = 0;
F.clauses{1}.cut = 0;
