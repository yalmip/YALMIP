function C = chanceconstraint(C,level)
C = flatten(C);
C.clauses{1}.jointprobabilistic = C.LMIid;
C.clauses{1}.confidencelevel = level;
for i = 2:length(C.clauses)
    C.clauses{i}.jointprobabilistic = C.LMIid;
    C.clauses{i}.confidencelevel = level;
end