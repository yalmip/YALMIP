function newConstraint = nonsymmetricUnivariateChanceFilter(b,c,distribution,gamma,w,options)

if isa(c,'sdpvar')
    error('Decision variables in c not yet supported')
else
    if length(c)==1
        if c>=0
            % Such as -t+w >= 0 i.e. t <= w
            newConstraint = -b <= c*icdf(distribution.parameters{1},gamma,distribution.parameters{2:end});
        else
            newConstraint = b >= -c*icdf(distribution.parameters{1},1-gamma,distribution.parameters{2:end});
        end
    else
        switch distribution.parameters{1}
            case 'gamma'
                k = distribution.parameters{2};               
                theta = distribution.parameters{3};               
                phi = @(t)(1-i*theta(:).*t).^-k;                               
                phi = @(t) prod(phi(-c(:)*t),1);
                newConstraint = [cdf_from_characteristic(b,phi) >= 1-gamma];
            otherwise
                error
        end
    end
end