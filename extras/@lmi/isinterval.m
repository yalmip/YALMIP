function YESNO = isinterval(F)
%ISINTERVAL (overloaded)

F = flatten(F);
if isempty(F.clauses)
    YESNO = 1;
else
    YESNO = 1;
    i = 1;
    while (i<=length(F.clauses)) & YESNO
        Fi = F.clauses{i};
        YESNO =  YESNO & is(Fi.data,'interval');
        i = i+1;
    end   
end

