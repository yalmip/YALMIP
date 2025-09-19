function sampledData = dataSampler(distribution,dimData)

if ~isempty(distribution.mixture)
    sampledData = mixLayer(distribution,dimData);
    return
end

temp = {distribution.generator,distribution.parameters{:},dimData};
if any(cellfun('isclass',temp,'sdpvar'))
    error('Cannot sample when parameters are decision variables')
end
if strcmp(func2str(temp{1}),'random')
    % We're most likely using a standard distribution from
    % the random command. However, we must check for some
    % additions to the standard cases
    if strcmp(temp{2},'normal')         
        % Normal is currently horrible. In some parts of the code, we
        % extend the parameter list to use both std and cov
        if length(distribution.parameters) == 4
            stD = distribution.parameters{end};
            if size(stD,2)>1
                stD = diag(stD);
            end
            temp = {distribution.generator,distribution.parameters{1:end-1},dimData};
            temp{4} = stD;         
        end
        sampledData = feval(temp{:});   
    elseif strcmp(temp{2},'mvnrnd')              
        sampledData = mvnrnd(temp{3}, temp{4});
        sampledData = sampledData(:);        
    elseif strcmp(temp{2},'mvnrndfactor')        
        sampledData = mvnrnd(temp{3}, temp{4}'*temp{4});
        sampledData = sampledData(:);        
    elseif strcmp(temp{2},'laplace')                
        sampledData1 =  random('exponential',temp{4}/sqrt(2));
        sampledData2 =  random('exponential',temp{4}/sqrt(2));
        sampledData = sampledData1 - sampledData2 + temp{3};
    elseif strcmp(temp{2},'data')
        i = randi(size(temp{3},2));
        sampledData = temp{3}(:,i);        
    else
        % sample from built-in standard distribution
        sampledData = feval(temp{:});   
    end
else
    sampledData = feval(temp{:});
end


function sampledData = mixLayer(distribution,dimData)

% Draw a mixture for every component
component = [];
for i = 1:size(distribution.mixture,1)
    n = size(distribution.mixture,2);
    cdf = cumsum(distribution.mixture(i,:));
    u = rand();
    component = [component;find(u <= cdf, 1)];
end
D = distribution;D.mixture = [];
for i = 2:length(distribution.parameters)
    p = [];
    for j = 1:size(distribution.parameters{i}{1},1)
        S = D.parameters{i}{component(j)};if size(S,2)>1;S = diag(S);end
        p = [p;S(j)];
    end
    D.parameters{i} = p;
end
sampledData = dataSampler(D,dimData);

