function newConstraint = symmetricUnivariateChanceFilter(b,c,distribution,gamma,options,funcs,x)

% General case must be handled via characteristic framework
if length(c) > 1 || trueishSetting(options.chance.characteristic)
    newConstraint = [(characteristic_cdf(x,funcs,distribution) >= 1-gamma) : 'Chance constraint (characteristics)'];
    return
end

% Scalar case, so we can derive an analytical condition
switch lower(distribution.parameters{1})
    
    case 't'                
        newConstraint = [b >= abs(c)*icdf(distribution.parameters{1},1-gamma,distribution.parameters{2:end})];
        
    case 'logistic'
        mu = distribution.parameters{2};
        distribution.parameters{2} = mu*0;
        newConstraint = b + c*mu >= abs(c)*icdf(distribution.parameters{1},1-gamma,distribution.parameters{2:end});
        
    case 'uniform'
        L = distribution.parameters{2};
        U = distribution.parameters{3};
        theMean = (L+U)/2;
        theRange = (U-L)/2;
        b = b+c*theMean;
        c = c*theRange;
        newConstraint = [b >= (1-2*gamma)*abs(c)];
        
    case 'laplace'
        mu = distribution.parameters{2};
        sigma = distribution.parameters{3};
        scale = sigma/sqrt(2);
        % Convex part of cdf for interesting high fidelity part simple!
        % b = (b+c'*mu)/abs(c);
        % cdf = @(x)(0.5*(1 + (1 - exp(-x/scale))));
        % newConstraint = [cdf(b) >= 1-gamma, 1 >= 1-gamma >= 1/2];
        newConstraint = [b+c*mu + abs(c)*scale*log(2*gamma) >= 0, 1 >= 1-gamma >= 1/2];
        
  case 'cauchy'
        mu = distribution.parameters{2};
        theta = distribution.parameters{3};        
        newConstraint = [b + c*mu >= abs(c)*theta*tan(pi*((1-gamma)-1/2))];        
        
    otherwise
        % Assume first parameter is mean (it was sent here, so it is
        % symmetric and thus first arument should be mean, right?...)
        mu = distribution.parameters{2};
        distribution.parameters{2} = mu*0;
        newConstraint = b + mu*c >= abs(c)*icdf(distribution.parameters{1},1-gamma,distribution.parameters{2:end});
        
end
newConstraint = [newConstraint : 'Chance constraint (analytic)'];