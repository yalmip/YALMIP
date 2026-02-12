function newConstraint = normalChanceFilter(b,c,distribution,gamma,options,funcs,x,isDisjointProblem)

% Formulating P(b(x) + c(x)'*w >= 0) >= 1-gamma w with Gaussian 

if length(b) > 1
    
    % If user inserted constraints without random variable, remove columns
    % of the corresponding constraint in c
    if isa(c,'sdpvar')
        % Not used at the moment but might be useful in the future
        baseC = getbase(c);
        colsAllZero = all(baseC==0,1);
    else
        colsAllZero = all(c==0,1);
    end
    if any(colsAllZero)
        c(:,colsAllZero) = [];
        b(colsAllZero) = [];
    end

    % The covariance
    Sigma = distribution.parameters{4};
    SigmaX = c'*Sigma*c;
    SigmaX = (SigmaX + SigmaX')/2; % Ensure numerical symmetry

    % Eigenvalue decomposition
    [theta,eigValues] = eig(SigmaX); % theta'*SigmaX*theta = eigvalues
    lambda = diag(eigValues);
    lambda(lambda<0) = 0; % Ensure lambdas are positive (numerical noise might give -1e-16)
    nphi = length(lambda); % Dimension of random vector phi
    nm = size(SigmaX,2); % Number of constraints

    % check orthogonality and diagonalization
    assert(norm(theta'*theta-eye(size(theta,2)),'inf') < 1e-12,'Theta not orthonormal within tol 1e-8');
    assert(norm(theta'*SigmaX*theta-eigValues,'fro') < 1e-12,'Diagonalization failed within tol 1e-8');

    
    % newConstraint = normalChanceFilterJointLayer(b,c,distribution,gamma,options,funcs,x,isDisjointProblem);
    % return
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
% Trivial Bonferroni
m = length(b);
gamma_local = sdpvar(m,1);
newConstraint = sum(gamma_local) <= gamma;
for i = 1:m
    newConstraint = [newConstraint, normalChanceFilter(b(i),c(:,i),distribution,gamma_local,options,funcs,x,isDisjointProblem)];
end
   
function C = allSequences(m,n)
% Returns an (m^n) x n matrix whose rows are all length-n sequences from 1:m.
vectors = cell(1,n);
[vectors{:}] = ndgrid(1:m);
C = cell2mat( cellfun(@(V) V(:), vectors, 'UniformOutput', false) );

