function [F,changed] = convertsocp(F)
%catsdp          Internal function: converts SOCP to LMI

% Author Johan Löfberg 
% $Id: convertsocp.m,v 1.3 2005-02-04 10:10:26 johanl Exp $

changed = 0;
Counter = size(F.clauses,2);
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

