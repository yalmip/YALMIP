function id = cleardual(F)

for i = 1:length(F.LMIid)
    yalmip('cleardual',F.LMIid(i));
end
