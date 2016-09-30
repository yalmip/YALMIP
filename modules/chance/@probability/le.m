function ChanceConstraint = le(level,P)

if isa(P,'double') & isa(level,'probability')
    error('Currently only supports level <= p(F)')
    temp = P;
    P = level;
    level = temp;
end

if level < 0 | level > 1
    error('The confidence level must be between 0 and 1');
end
ChanceConstraint = chance(P.Constraint,level);

