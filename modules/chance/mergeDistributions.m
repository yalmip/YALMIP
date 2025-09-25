function [merged,map] = mergeDistributions(randomVariables)

merged = {};
map = [];
used = zeros(1,length(randomVariables));
mergedindex = 1;

% Mainly normalizing Gaussian variables as they come in 3 different ways
for i = 1:length(randomVariables)        
	randomVariables{i} = normalize(randomVariables{i});
end

for i = 1:length(randomVariables)
    if ~used(i)
        this = randomVariables{i};        
        used(i) = 1;
        for j = i+1:length(randomVariables)
            if ~used(j)
                that = randomVariables{j};
                if strcmp(func2str(this.distribution.generator),'random') && strcmp(func2str(that.distribution.generator),'random')
                    if mergable(this.distribution.parameters{1},that.distribution.parameters{1})
                        % Same distribution                        
                        this = merge2distributions(this,that);
                        map(j) = mergedindex;
                        used(j) = 1;
                    end
                end
            end
        end
        merged{mergedindex} = this;
        map(i) = mergedindex;
        mergedindex = mergedindex + 1;
    end
end


function normalized = normalize(D)

% A bit messy as univariate and multivariate are defined differently in
% statistics toolbox and then we have a YALMIP specific also
% normal, scalar elementwise normal (parameterized in std. dev!)
% mvrnd, multivariable normal, (parameterized in covariance!)
% mvrndfactor, m.v. normal, (parameterized in factor covariance)

% To simplify merging, we introduce generalized 'normal' and equip it
% with both covariance and factorized covariance

normalvariants = {'normal','mvnrnd','mvnrndfactor'};
normalized = D;
if strcmp(func2str(D.distribution.generator),'random')
    if any(strcmp(D.distribution.parameters{1},normalvariants))                       
        if isempty(D.distribution.mixture)
            [varD,facD] = normalizeOneTerm(D.distribution.parameters{1},D.distribution.parameters{3},D.variables);
            normalized = D;
            normalized.distribution.parameters{1} = 'normal';
            normalized.distribution.parameters{3} = varD;
            normalized.distribution.parameters{4} = facD;         
        else
            for i = 1:size(D.distribution.mixture,2)
                [varD{i},facD{i}] = normalizeOneTerm(D.distribution.parameters{1},D.distribution.parameters{3}{i},D.variables);
            end
            normalized.distribution.parameters{1} = 'normal';
            normalized.distribution.parameters{3} = varD;
            normalized.distribution.parameters{4} = facD;
        end
    end
end

function [varD,facD] = normalizeOneTerm(dName,spreadParam,variables)

dimDvar = size(spreadParam);
if dimDvar(1) ~= dimDvar(2)
    parD = diag(spreadParam);
else
    parD = spreadParam;
    if numel(parD)==1 && numel(variables)>1
        parD = diag(repmat(parD,numel(variables),1));
    end
end
if strcmp(dName,'normal') || strcmp(dName,'mvnrndfactor')
    varD = parD'*parD;
    facD = parD;
elseif strcmp(dName,'mvnrnd')
    varD = parD;
    if isa(parD,'double')
        facD = chol(parD);
    else
        facD = [];
    end
end
