function LinearVariables = depends(F)

% Get all used variables in this LMI object
used = [];
for i = 1:length(F.clauses)  
    Fivar = getvariables(F.clauses{i}.data);
    used = uniquestripped([used Fivar(:)']);
end

% Now, determine the involved linear variables
[mt,variabletype] = yalmip('monomtable');
if any(variabletype)%1%~isempty(nlv) & any(ismembc(used,nlv(1,:)))
    LinearVariables = find(any(mt(used,:),1));    
    LinearVariables = LinearVariables(:)';
else
    LinearVariables = used;
end
