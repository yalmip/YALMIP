function sampledData = dataSampler(distribution,dimData)
temp = {distribution.generator,distribution.parameters{:},dimData};
if any(cellfun('isclass',temp,'sdpvar'))
    error('Cannot sample when parameters are decision variables')
end
if strcmp(func2str(temp{1}),'random')
    % We're most likely using a standard distribution from
    % the random command. However, we must check for some
    % additions to the standard cases
    if strcmp(temp{2},'mvnrnd')              
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