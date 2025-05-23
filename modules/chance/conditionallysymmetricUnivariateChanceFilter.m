function newConstraint = conditionallysymmetricUnivariateChanceFilter(b,c,distribution,gamma,w,options,funcs,x)

% General case must be handled via characteristic framework
if length(c) > 1 || strcmpi(options.chance.characteristic,'yes')
    newConstraint = [characteristic_cdf(x,funcs,distribution) >= 1-gamma];
    return
end

% Scalar case, so we can derive an analytical condition
switch lower(distribution.parameters{1})
    
    case 'stable'                
        if distribution.parameters{3} == 0
            mu = distribution.parameters{5};
            distribution.parameters{5} = 0;
            newConstraint = b + c'*mu >= abs(c)*icdf(distribution.parameters{1},1-gamma,distribution.parameters{2:end});
        else
            newConstraint = [characteristic_cdf(x,funcs,distribution) >= 1-gamma];
        end
        
    otherwise
       error('')
        
end