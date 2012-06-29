function ChanceConstraint = gt(P,level)
if isa(P,'double') & isa(level,'probability')
    error('Currently only supports p(F) >= level')
    %temp = P;
    %P = level;
    %level = temp;
end
if level < 0 | level > 1
    error('The confidence level must be between 0 and 1');
end
ChanceConstraint = chanceconstraint(set(P.Constraint),level);

