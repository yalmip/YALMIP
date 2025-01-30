function newConstraint = symmetricUnivariateChanceFilter(b,c,distribution,gamma,w,options,funcs,x)

switch lower(distribution.parameters{1})
    
    % We have some special-case code for some distributions
    case 'uniform'              
        if length(c)==1
            % icdf is trivial
            % Convert to symmetric around origin for simplicity
            L = distribution.parameters{2};
            U = distribution.parameters{3};            
            theMean = (L+U)/2;
            theRange = (U-L)/2;
            b = b+c'*theMean;
            c = c.*theRange;
            newConstraint = [b >= (1-2*(1-gamma))*abs(c)];
        else
            phi = distribution.characteristicfunction;            
            dphi = distribution.characteristicfunction_derivative;                                    
            newConstraint = [characteristic_cdf(x,funcs,phi,dphi) >= 1-gamma];           
        end
        
    case 'laplace'
        mu = distribution.parameters{2};
        sigma = distribution.parameters{3};
        % Alternative notion with scale parameter common
        scale = sigma/sqrt(2);
        % We're shfting mean to 0 so numerator becomes 1 in characteristic
        phi = @(t) 1./(1+scale(:).^2.*t.^2);
        
        if length(c)==1 && strcmpi(options.chance.expcone,'yes')
            % Convex part of cdf for interesting high fidelity part simple!
            b = (b+c'*mu)/c;
            cdf = @(x)(0.5*(1 + (1 - exp(-x/scale))));
            newConstraint = [cdf(b) >= 1-gamma, 1 >= 1-gamma >= 1/2];
        elseif strcmpi(options.chance.expcone,'yes')
            error('Exponential cone representation not available')
        else
            b = (b+c'*mu);
            phi = @(t) prod(phi(c(:)*t),1);
            newConstraint = [cdf_from_characteristic(b,phi) >= 1-gamma];
        end
        
    case 'logistic'
        mu = distribution.parameters{2};
        s = distribution.parameters{3};
        if length(c) == 1
            newConstraint = b >= abs(c)*icdf(distribution.parameters{1},1-gamma,distribution.parameters{2:end});
            % Possible to exploit on region of interest?
            % newConstraint = b   >= abs(c)*(mu+s*(log((1-gamma))-log(gamma)));
        else            
            phi = distribution.characteristicfunction;            
            dphi = distribution.characteristicfunction_derivative;            
            
            % Testing
            newConstraint = [characteristic_cdf(x,funcs,phi,dphi) >= 1-gamma];                        
        end
        
    otherwise
        newConstraint = b >= icdf(distribution.parameters{1},1-gamma,distribution.parameters{2:end})*abs(c);
end