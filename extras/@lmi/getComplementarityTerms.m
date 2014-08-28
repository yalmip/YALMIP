function [x,y] = getComplementarityTerms(F)

F = flatten(F);
xy = F.clauses{1}.data;
x = xy(:,1);
y = xy(:,2);
