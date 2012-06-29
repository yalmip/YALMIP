function [F,changed] = convertlorentz(F)
%catsdp          Internal function: converts rotated Lorentz to SOCC
  
% Author Johan Löfberg 
% $Id: convertlorentz.m,v 1.4 2005-09-23 13:05:00 joloef Exp $
  
changed = 0;
Counter = size(F.clauses,2);
for i = 1:Counter
    % Yep, Lorentz
    if  (F.clauses{i}.type==5)
        changed=1;
        xyz = F.clauses{i}.data;
        x = xyz(1);
        y = xyz(2);
        z = xyz(3:end);
        F.clauses{i}.data = [(x+y)/sqrt(2);(x-y)/sqrt(2);z];
        F.clauses{i}.type = 4;
    end
end

