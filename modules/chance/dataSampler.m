function sampledData = dataSampler(distribution,dimData);
temp = {distribution.name,distribution.parameters{:},dimData};
if strcmp(func2str(temp{1}),'random')
    % We're most likely using a standard distribution from
    % the random command. However, we must check for some
    % additions to the standard cases
    % 1. 'normalf' distribution, factorized matrix covariance
    % 2. 'normalm' distribution, matrix covariance
    % 3. sdpvar objects in arguments
    if strcmp(temp{2},'normal')
        sampledData = feval(temp{1:2},0,1,temp{5});
        sampledData = sqrtm(temp{4})*sampledData;
        sampledData = sampledData + temp{3};
    elseif strcmp(temp{2},'normalf')
        temp{2} = 'normal';
        sampledData = feval(temp{1:2},0,1,temp{5});
        sampledData = temp{4}*sampledData;
        sampledData = sampledData + temp{3};
    elseif strcmp(temp{2},'normalm')
        temp{2} = 'normal';
        sampledData = feval(temp{1:2},0,1,temp{5});
        sampledData = chol(temp{4})*sampledData;
        sampledData = sampledData + temp{3};
    elseif strcmp(temp{2},'data')
        i = randi(size(temp{3},2));
        sampledData = temp{3}(:,i);
    elseif ~any(cellfun('isclass',temp,'sdpvar'))
        % simple numerical sample of standard distriution
        sampledData = feval(temp{:});
    else
        error('Cannot sample from this description')
    end
else
sampledData = feval(temp{:});
end