function self = sample(self,N)
%SAMPLE Draw a sample in an optimizer object and instantiates parameter
%
%   Q = SAMPLE(P,N) generates concatenated instantiated optimizer objects
%   where samples have been drawn for all parameters with associated
%   samplers.
%
%   INPUT
%    P : OPTIMIZER object
%    N : Number of samples to draw
%
%   OUTPUT
%    Q : OPTIMIZER object
%
%   EXAMPLE
%    sdpvar x(2,1) w(2,1)
%    F = [-10 <= x + w <= 10, uncertain(w,'unif',2,3)]
%    P = optimizer(F,sum(x),[],w,x)
%    Q = sample(P,10);
%    plot(Q);
%    xoptimal = Q([]);
%
%   See also OPTIMIZER, UNCERTAIN

if nargin < 2
    N = 1;
end
allSamples = {};
for k = 1:N
    cells = cell(1,length(self.diminOrig));
    for i = 1:length(self.diminOrig)
        if isfield(self.input,'stochastics')
            if ~isempty(self.input.stochastics{i});
                temp = {self.input.stochastics{i}.name,self.input.stochastics{i}.parameters{:},self.diminOrig{i}};
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
                        sampledData = feval(temp{1:2},0,1,temp{5});
                        sampledData = temp{4}*sampledData;
                        sampledData = sampledData + temp{3};
                    elseif strcmp(temp{2},'normalm')
                        sampledData = feval(temp{1:2},0,1,temp{5});
                        sampledData = chol(temp{4})*sampledData;
                        sampledData = sampledData + temp{3};
                    elseif ~any(cellfun('isclass',temp,'sdpvar'))
                        % simple numerical sample of standard distriution
                        sampledData = feval(temp{:});
                    else
                        error('Cannot sample from this description')
                    end
                else
                    sampledData = feval(temp{:});
                end
                cells{i} = sampledData;
            else
                cells{i} = [];
            end
        end
    end
    if ~isempty(cells) && ~all(cellfun('isempty',cells))
        y.type = '{}';
        y.subs = {cells{:},'nosolve'};      
        allSamples{end + 1} = subsref(self,y);
    else
        warning('Trying to sample in a model without any random uncertainties');
        return
    end
end
self = horzcat(allSamples{:});


