function lrdata = getlrdata(F)

k = 1;
lrdata = [];
for i = 1:length(F.clauses)
    if F.clauses{i}.type == 14
        lrdata{end+1} =  F.clauses{i}.data;
    end    
end
