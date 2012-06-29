function B = getbase(F)

if length(F.clauses)>1
    error('LMI/GETBASE can only be applied to a list with 1 constraint')
else
    B = getbase(F.clauses{1}.data);
end
