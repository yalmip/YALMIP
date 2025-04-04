function ChanceConstraint = le(level,P)

if isa(P,'double') & isa(level,'probability')
    error('Currently only supports [probability(f) >= level]')  
end

if isa(level,'double') && level < 0 || level > 1
    error('The confidence level must be between 0 and 1');
end

ChanceConstraint = chanceconstraint(lmi(P.Constraint),level);

