function YESNO = isinterval(F)
%ISINTERVAL (overloaded)

% Author Johan Löfberg 
% $Id: isinterval.m,v 1.1 2007-03-23 12:36:02 joloef Exp $   

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

