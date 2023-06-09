function [F] = registerjacobian(F)
%REGISTERJACOBIAN Register jacobians to all constraints

F = flatten(F);
Counter = length(F.LMIid);
for i = 1:Counter
    switch F.clauses{i}.type
        case {1,2,3}
            Fi = F.clauses{i}.data;
            F.clauses{i}.data = registerjacobian(Fi);
        otherwise
    end
end
