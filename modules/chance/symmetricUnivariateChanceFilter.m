function newConstraint = symmetricUnivariateChanceFilter(b,c,distribution,gamma,w,options)

switch distribution.parameters{1}
    
    case 'uniform'
        % not necessarily symmetric, but we can shift it
        L = distribution.parameters{2};
        U= distribution.parameters{3};
        theMean = (L+U)/2;
        theRange = (U-L)/2;
        b = b+c*theMean;
        c = c*theRange;
        distribution.parameters{2}=-1;
        distribution.parameters{3}= 1;
        
        % and icdf is trivial
        if isa(gamma,'sdpvar')
            newConstraint = [b >= (1-2*gamma)*abs(c), 1>=gamma>=0];
        else
            newConstraint = [b >= (1-2*gamma)*abs(c)];
        end
        
    otherwise
        newConstraint = b >= icdf(distribution.parameters{1},1-gamma,distribution.parameters{2:end})*abs(c);
end