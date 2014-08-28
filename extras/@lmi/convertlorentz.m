function [F,changed] = convertlorentz(F)
%convertlorentz          Internal function: converts rotated Lorentz to SOCC
  
changed = 0;
F = flatten(F);
Counter = length(F.LMIid);
for i = 1:Counter
    % Yep, Lorentz
    if  (F.clauses{i}.type==5)
        changed=1;
        xyz = sdpvar(F.clauses{i}.data);
        x = xyz(1);
        y = xyz(2);
        z = xyz(3:end);
        F.clauses{i}.data = [(x+y)/sqrt(2);(x-y)/sqrt(2);z];
        F.clauses{i}.type = 4;
    end
end

