function F = remap(F,old,new)

for i = 1:length(F.clauses)    
    X = F.clauses{i}.data;    
    X = sdpvarremap(X,old,new);    
    F.clauses{i}.data = X;
end