function data = sample(w,N)

if nargin == 1
    N = 1;
end
randomVariables = yalmip('getDistribution');
% Prune unused random declarations to avoid extra work. ALso check for
% special case where w is a simple variable
allVars = depends(w);
allwVars = [];
keep = false(1,length(randomVariables));
simpleCase = [];
for i = 1:length(randomVariables)
    wVars = getvariables(randomVariables{i}.variables);
    if any(ismember(allVars, wVars))
        allwVars = [allwVars;wVars];
        keep(i) = true;
        if isequal(getvariables(w),getvariables(randomVariables{i}.variables))
            if isequal(getbase(w), getbase(randomVariables{i}.variables))
                simpleCase = i;
                break
            end
        end
    end
end

if ~simpleCase
    randomVariables = randomVariables(keep);
end

if length(randomVariables) == 0
    data = repmat(w,1,N);
    return
end

if ~isempty(simpleCase)
    data = [];
    distribution = randomVariables{simpleCase}.distribution;
    data = getSamples(distribution,N);
else
    error
end

function data = getSamples(distribution,N)

switch func2str(distribution.generator)
    case 'random'
        switch distribution.parameters{1}
            case 'mvnrnd'
                data = [];
                for i = 1:N
                    s = mvnrnd(distribution.parameters{2:end});
                    data = [data s(:)];
                end
            case 'data'
                i = randi(size(distribution.parameters{2},2),N,1);
                data = distribution.parameters{2}(:,i);
            otherwise
                data = [];
                for i = 1:N
                    data = [data random(distribution.parameters{:})];
                end
        end
    otherwise
        error
end