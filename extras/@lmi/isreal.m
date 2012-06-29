function YESNO = isreal(F)
%ISREAL (overloaded)

% Author Johan Löfberg 
% $Id: isreal.m,v 1.3 2007-03-23 12:36:02 joloef Exp $   

if isempty(F.clauses)
    YESNO = 1;
else
    YESNO = 1;
    i = 1;
    while (i<=length(F.clauses)) & YESNO
        Fi = F.clauses{i};
        YESNO =  YESNO & is(Fi.data,'real');
        i = i+1;
    end   
end

