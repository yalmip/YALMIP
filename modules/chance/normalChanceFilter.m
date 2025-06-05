function newConstraint = normalChanceFilter(b,c,distribution,gamma,options,funcs,x,isDisjointProblem)

% Mixture layer
if ~isempty(distribution.mixture)
    newConstraint = normalChanceFilterMixLayer(b,c,distribution,gamma,options,funcs,x,isDisjointProblem);
    return
end

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
        error('Cannot have decision variables multplying uncertainty when using expcone approximation of inverse cdf')
    end
    Phi_Inverse = normalChanceFilterConicApproximation(gamma);
    newConstraint =  b + c'*theMean >= Phi_Inverse*norm(e);
elseif strcmpi(options.chance.expcone,'root') && isDisjointProblem
    rootPhi_Inverse = normalChanceFilterConicApproximationRoot(gamma);
    newConstraint = normalChanceFilterConicFormulationRoot(c,b,rootPhi_Inverse);
elseif strcmpi(options.chance.expcone,'log') && isDisjointProblem
    logPhi_Inverse = normalChanceFilterConicApproximationLog(gamma);
    newConstraint = normalChanceFilterConicFormulationLog(c,b,logPhi_Inverse);
elseif strcmpi(options.chance.expcone,'inv') && isDisjointProblem
    invPhi_Inverse = normalChanceFilterConicApproximationInv(gamma);
    newConstraint = normalChanceFilterConicFormulationInv(c,b,invPhi_Inverse);
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



function newConstraint = normalChanceFilterMixLayer(b,c,distribution,gamma,options,funcs,x,isDisjointProblem)
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
factorcovariance = distribution.parameters{4};

% FIXME: Refactor
% Mixtures are currently only acting on elementwise terms. This means
% that the mixture model now has to be exploded to be the mixture of
% the linear combination
newWeights = [];
newMean = {};
newcovariance = {};
n = length(c);
m = length(distribution.mixture);
Combs = allSequences(m,n);
newMixture = [];
for i = 1:length(Combs)
    if any(cellfun(@(x)isa(x,'sdpvar'),theMean))
        aMean = sdpvar(n,1);
    else
        aMean = zeros(n,1);
    end
    aCov =  zeros(n,1);
    for j = 1:n
        aMean(j) = theMean{Combs(i,j)}(j);
        aCov(j) = covariance{Combs(i,j)}(j,j);
    end
    allMeans{i} = aMean;
    allCovs{i} = diag(aCov);
    newMixture(end+1) = prod(distribution.mixture(Combs(i,:)));
end
distribution.mixture = [];
gamma_i = sdpvar(1,length(newMixture));
newConstraint = [gamma == sum(newMixture.*gamma_i)];
for i = 1:length(newMixture)
    distribution.parameters{2} = allMeans{i};
    distribution.parameters{3} = allCovs{i};
    distribution.parameters{4} = allCovs{i}.^.5;
    newConstraint = [newConstraint, normalChanceFilter(b,c,distribution,gamma_i(i),options,funcs,x,isDisjointProblem)];
end
return



function C = allSequences(m,n)
% Returns an (m^n) x n matrix whose rows are all length-n sequences from 1:m.
vectors = cell(1,n);
[vectors{:}] = ndgrid(1:m);
C = cell2mat( cellfun(@(V) V(:), vectors, 'UniformOutput', false) );

