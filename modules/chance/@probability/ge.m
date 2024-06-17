function ChanceConstraint = gt(P,level)

if isa(P,'double') & isa(level,'probability')
    error('Currently only supports [probability(f) >= level')
end

if isa(level,'double') && (level < 0 || level > 1)
    error('The confidence level must be between 0 and 1');
end

if length(P.Constraint) > 1
    % Rewrite P(C1) + P(C2) + ... + P(Cn) >= 1-gamma = level as
    % Rewrite P(Ci) >= 1-risk_i, sum 1-risk_i >= level
    ChanceConstraint = [];
    lhs = 0;
    for i = 1:length(P.Constraint)
        if isempty(P.Constraint{i})
            lhs = lhs + P.Offset{i};
        else
            lhs = lhs + (1-P.Risk{i});
            ChanceConstraint = [ChanceConstraint,chanceconstraint(lmi(P.Constraint{i}),1-P.Risk{i})];
        end
    end
    ChanceConstraint = [ChanceConstraint, lhs >= level];
else
    ChanceConstraint = chanceconstraint(lmi(P.Constraint{1}),level);
end

