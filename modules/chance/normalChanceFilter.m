function newConstraint = normalChanceFilter(b,c,distribution,gamma,options,funcs,x,isDisjointProblem)

% Formulating P(b(x) + c(x)'*w >= 0) >= 1-gamma w with Gaussian

if length(b) > 1

    newConstraint = normalChanceFilterJointLayer(b,c,distribution,gamma,options,funcs,x,isDisjointProblem);
    return
end

if strcmpi(options.chance.characteristic,'yes')
    % A bit messy with normal/gaussian. Internal framework works with mean,
    % covariance and factorized covariance (std. dev in scalar case) and
    % thus has 3 parameters. However, the characteristic function is
    % defined from mean and std. dev to comply with standard notation.
    % Hence remove second parameter and extract diagonal of factor
    distribution.parameters = {distribution.parameters{1:2}, cellfun(@(c)diag(c),distribution.parameters{4},'UniformOutput', false)};
    % distribution.parameters = {distribution.parameters{1:2}, distribution.parameters{4}};
    newConstraint = [characteristic_cdf(x,funcs,distribution) >= 1-gamma];
    return
end

% default though is to simply use analytic stuff
if ~isempty(distribution.mixture)
    newConstraint = normalChanceFilterMixLayer(b,c,distribution,gamma,options,funcs,x,isDisjointProblem);
    return
end

theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
factorcovariance = distribution.parameters{4};

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



function newConstraint = normalChanceFilterMixLayer(b,c,distribution,gamma,options,funcs,x,isDisjointProblem)
theMean    = distribution.parameters{2};
covariance = distribution.parameters{3};
factorcovariance = distribution.parameters{4};

% There might be zeros in c meaning some w(i) are missing and can be pruned
missing = ~any(getbase(c),2);
if ~isempty(missing)
    c(missing)=[];
    for i = 1:length(theMean)
        theMean{i}(missing) = [];
        covariance{i}(:,missing) = [];
        covariance{i}(missing,:) = [];
        factorcovariance{i}(:,missing) = [];
    end
end
% FIXME: Refactor
% Mixtures are currently only acting on elementwise terms. This means
% that the mixture model now has to be exploded to be the mixture of
% the linear combination
newWeights = [];
newMean = {};
newcovariance = {};
n = length(c);
m = size(distribution.mixture,2);
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
    mix = 1;
    for k = 1:size(Combs,2);mix = mix*distribution.mixture(k,Combs(i,k));end
    newMixture(end+1) = mix;
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


function newConstraint = normalChanceFilterJointLayer(b,c,distribution,gamma,options,funcs,x,isDisjointProblem)

% The covariance
factorcovariance = distribution.parameters{4};
Sigma = c'*factorcovariance*c;
Sigma = (Sigma + Sigma')/2; % Ensure numerical symmetry

% Eigenvalue decomposition
[theta,eigValues] = eig(Sigma); % theta'*SigmaX*theta = eigvalues
lambda = diag(eigValues);
lambda(lambda<0) = 0; % Ensure lambdas are positive (numerical noise might give -1e-16)

% check orthogonality and diagonalization
assert(norm(theta'*theta-eye(size(theta,2)),'inf') < 1e-12,'Theta not orthonormal within tol 1e-8');
assert(norm(theta'*Sigma*theta-eigValues,'fro') < 1e-12,'Diagonalization failed within tol 1e-8');

% Auxiliary betas
beta1 = sdpvar(length(lambda),1);
beta2 = sdpvar(length(lambda),1);

% Create the product constraint with log for computational burden
logBeta = 0;
for i = 1:length(lambda)
    logBeta = logBeta + log(beta1(i) + beta2(i) - 1);
end
newConstraint = [logBeta >= log(1-gamma),
    beta1 + beta2 >= 1,
    0 <= beta1 <= 1,
    0 <= beta2 <= 1];

for i = 1:size(Sigma,2)
    lhs_sum = 0;
    for j = 1:length(lambda)
        % Check sign of theta_ij to select beta1 or beta2
        if theta(i,j) >= 0
            b_val = beta1(j);
        else
            b_val = beta2(j);
        end
        lhs_sum = lhs_sum + sqrt(lambda(j))*abs(theta(i,j))*norminv(b_val);
    end

    newConstraint = [newConstraint, lhs_sum - b(i) <= 0];

end
return

