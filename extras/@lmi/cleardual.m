function id = cleardual(F)

for i = 1:length(F.clauses)
    yalmip('cleardual',F.LMIid(i));
end
