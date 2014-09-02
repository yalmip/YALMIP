function [F,changed] = convertsocp(F)
%catsdp          Internal function: converts SOCP to LMI

changed = 0;
F = flatten(F);
Counter = length(F.LMIid);
for i = 1:Counter
    if  (F.clauses{i}.type==4)
        changed=1;
        xy = F.clauses{i}.data;
        y = xy(1);
        x = xy(2:end);
        F.clauses{i}.data = [y x';x eye(length(x))*y];
        F.clauses{i}.type = 1;
    end
end

