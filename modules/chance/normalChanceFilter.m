function newConstraint = normalChanceFilter(b,c,distribution,gamma,w,options,funcs,x,isDisjointProblem)
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
factorcovariance = distribution.parameters{4};

% You never know, someone might want to do this...
if strcmpi(options.chance.characteristic,'yes') && isdiag(covariance)
    if isa(theMean,'double') && isa(covariance,'double')
        newConstraint = [characteristic_cdf(x,funcs,distribution) >= 1-gamma];
        return
    else
        error('Characterstics can only be used for fixed distribution parameters');
    end
end

constant_gain = isa(c,'double') && isa(factorcovariance,'double');

if constant_gain && (~strcmpi(options.chance.expcone,'no'))
    % c is constant, so no reason really to use log tricks etc?
    % for now, just revert to standard expcone approximation
    options.chance.expcone = 'yes';
end

e = factorcovariance*c;

if strcmpi(options.chance.expcone,'yes')
    if ~constant_gain
        error('Cannot have decision variables multiplying uncertainty when using expcone approximation of inverse cdf')
    end
    newConstraint = normalChanceFilterConicFormulation(b,c,e,theMean,gamma);
elseif strcmpi(options.chance.expcone,'inv') && isDisjointProblem
    newConstraint = normalChanceFilterConicFormulationInv(c,b,gamma);
elseif strcmpi(options.chance.expcone,'root') && isDisjointProblem
    newConstraint = normalChanceFilterConicFormulationRoot(c,b,gamma);
elseif strcmpi(options.chance.expcone,'log') && isDisjointProblem
    newConstraint = normalChanceFilterConicFormulationLog(c,b,gamma);
else
    % Just go for a general nonlinear model and hope for the best
    Phi_Inverse = icdf('normal',1-gamma,0,1);
    if isa(Phi_Inverse,'sdpvar')
        % If we have a general nonlinear model, we should not use norm as it
        % is intended to be socp-represented. Avoid that by using the callback
        % version of norm
        % Note, this cannot happen when expcone is used, stopped above
        newConstraint = b + c'*theMean >= Phi_Inverse*norm_callback(e);
    else
        newConstraint =  b + c'*theMean >= Phi_Inverse*norm(e);
    end
end