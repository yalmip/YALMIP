function ChanceConstraint = ge(P,level)

if isa(P,'double') & isa(level,'probability')
    error('Currently only supports [probability(f) >= level')
end

if isa(level,'double') && (any(level < 0) || any(level > 1))
    error('The confidence level must be between 0 and 1');
end

if length(P.Constraint) > 1
    % Rewrite P(C1) + P(C2) + ... + P(Cn) >= 1-gamma = level as
    % Rewrite P(Ci) >= 1-risk_i, sum 1-risk_i >= level
    ChanceConstraint = [];
    lhs = 0;
    for i = 1:length(P.Constraint)
        if isempty(P.Constraint{i})
            lhs = lhs + P.Weight{i}*P.Offset{i};
        else
            lhs = lhs + P.Weight{i}*(1-P.Risk{i});
            ChanceConstraint = [ChanceConstraint,chanceconstraint(lmi(P.Constraint{i}),1-P.Risk{i})];
        end
        if isa(P.Weight{i},'sdpvar')
            ChanceConstraint = [ChanceConstraint, P.Weight{i} >= 0];
        end
    end
    ChanceConstraint = [ChanceConstraint, lhs >= level];
else
    % We add a redundant constraint to make sure the risk-variable
    % associated here is available also is trivial cases, for the value
    % operator to be able to return something
    % FIXME: Appears to trigger numerical issues sometimes in fmincon
    %ChanceConstraint = [1-P.Risk{1} >= level/P.Weight{1}, chanceconstraint(lmi(P.Constraint{1}),level/P.Weight{1})];
    c = sdpvar(P.Constraint{1});
    if length(c)>1 && length(level)>1 &&  length(c) ~=  length(level)
        error('Dimension mismatch in probabilistic constraint')
    elseif length(c)>1 && length(level)==1
        level = repmat(level,length(c),1);
    elseif length(c)==1 && length(level)>1
        c = repmat(c,length(level),1);
    end
    if length(c) == 1
        ChanceConstraint = [chanceconstraint(lmi(P.Constraint{1}),level/P.Weight{1})];
    else
        ChanceConstraint = [];
        for i = 1:length(c)
            ChanceConstraint = [ChanceConstraint, chanceconstraint(lmi(c(i)),level(i)/P.Weight{1})];
        end
    end            
end

