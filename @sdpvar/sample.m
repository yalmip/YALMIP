function wrealization = sample(Expression,N)

if nargin < 2
    N = 1;
end

% Get all definitions
randomVariables = yalmip('getDistribution');
% Scalar definitions of same distribution bunched together
% This merges various normal/gaussian models, beware
randomVariables = mergeDistributions(randomVariables);

% Now prune for those that are actually asked for
w_variables = depends(Expression);
keep = zeros(1,length(randomVariables));
simple = 0;
for i = 1:length(randomVariables)
    wi = randomVariables{i}.variables;
    if any(ismember(w_variables,getvariables(wi)))
        keep(i) = 1;
        if isequal(w_variables,getvariables(wi))
            if isequal(getbase(Expression),getbase(wi))
                simple = i;
                break
            end
        end
    end
end
randomVariables = {randomVariables{find(keep)}};
W = [];
for i = 1:length(randomVariables)
    W = [W;randomVariables{i}.variables];
end

% Generate samples for all involved random variables
allSamples = [];
for i = 1:length(randomVariables)
    wi = randomVariables{i}.variables;
    samples{i} = [];
    for j = 1:N
        samples{i} =  [samples{i} dataSampler(randomVariables{i}.distribution,size(wi),randomVariables{i}.id)];
    end
    allSamples = [allSamples;samples{i}];
end

% And now map to the requested variable
if simple
    wrealization = allSamples;
else
    % FIXME: Do this faster
    wrealization = [];
    for i = 1:N
        wrealization = [wrealization replace(Expression,W,allSamples(:,i))];
    end
end
