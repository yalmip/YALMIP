function newConstraint = normalChanceFilter(b,c,distribution,gamma,options,funcs,x,isDisjointProblem)

if length(b) > 1
    newConstraint = normalChanceFilterJointLayer(b,c,distribution,gamma,options,funcs,x,isDisjointProblem);
    return
end

if any(strcmpi(options.chance.characteristic,{'yes','on'})) || isequal(options.chance.characteristic,1)
    % A bit messy with normal/gaussian. Internal framework works with mean,
    % covariance and factorized covariance (std. dev in scalar case) and
    % thus has 3 parameters. However, the characteristic function is
    % defined from mean and std. dev to comply with standard notation.
    % Hence remove second parameter and extract diagonal of factor
    if isa(distribution.parameters{4},'cell')
        % Mixture case
        distribution.parameters = {distribution.parameters{1:2}, cellfun(@(c)diag(c),distribution.parameters{4},'UniformOutput', false)};    
    else
        distribution.parameters = {distribution.parameters{1:2}, diag(distribution.parameters{4})};
    end  
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

