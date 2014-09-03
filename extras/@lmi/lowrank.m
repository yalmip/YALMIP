function F = lowrank(F,x)
% LOWRANK is used to declare that a semidefinite constraint uses data with low rank.
% Used in combination with the solver SDPLR

identifiers = [];
F = flatten(F);
for i = 1:length(F.LMIid)
    switch F.clauses{i}.type
        case 1
            identifiers = [identifiers F.LMIid(i)];
        otherwise            
    end
end

if isempty(identifiers)
    F = ([]);
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
