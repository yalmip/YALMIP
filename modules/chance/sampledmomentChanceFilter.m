function newConstraint =  sampledmomentChanceFilter(b,c,distribution,confidencelevel,w,options);
theMean = [];
theCov = [];
if strcmpi(func2str(distribution.name),'random') && strcmpi(distribution.parameters{1},'data')
    % Only data given, so best we can do is to estimate from this
    W = distribution.parameters{2};
    N = size(W,2);  
    theMean = mean(W,2);
    theCov = cov(W');
    %d.parameters{2} = mean(W,2);
    %d.parameters{3} = cov(W');
elseif strcmpi(func2str(distribution.name),'random') 
    % We can compute the required moments
    [theMean, theCov] = extractMoments(distribution);
    if length(w) > 1 && length(distribution.parameters{2})==1
        theMean = repmat(theMean,length(w),1);
        theCov = theCov*eye(length(w));
    end
   % d.parameters{2} = theMean;
   % d.parameters{3} = theCov;
end
if isempty(theMean) || isempty(theCov)
    % General case, so draw samples to estimate mean and covariance
    % TODO: Don't though away analytical if available
    N = options.chance.N;
    W = [];for i = 1:N;W = [W dataSampler(distribution,size(w))];end    
    theMean = mean(W,2);
    theCov = cov(W');
end
d.parameters{2} = theMean;
d.parameters{3} = theCov;
newConstraint = momentChanceFilter(b,c,d,confidencelevel,w,options);

function [theMean, theCov] = extractMoments(distribution)
switch lower(distribution.parameters{1})
    case 'uniform'
        theMean = (distribution.parameters{3} + distribution.parameters{2})/2;
        theCov = (distribution.parameters{3} - distribution.parameters{2}).^2/12;
    case 'exponential'
        theMean = distribution.parameters{2};
        theCov = distribution.parameters{2}^2;
    otherwise
        theMean = [];
        theCov = [];
end