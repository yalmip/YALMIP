function newConstraint = nonsymmetricUnivariateChanceFilter(b,c,distribution,gamma,w,options)

if isa(c,'sdpvar')
    error('Decision variables in c not yet supported')
else
    if c>=0
        newConstraint = b >= -icdf(distribution.parameters{1},gamma,distribution.parameters{2:end})*c;
    else
        newConstraint = b >= -icdf(distribution.parameters{1},1-gamma,distribution.parameters{2:end})*c;
    end
end