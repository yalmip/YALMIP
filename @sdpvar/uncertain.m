function x = uncertain(x,varargin)
%UNCERTAIN Declares a variable as uncertain
%
%   F = UNCERTAIN(W) declares W to be (deterministic) uncertainty
%
%   F = UNCERTAIN(W,'name', param1, param2,...) declares W to be
%   stochastic uncertainty of type 'name' with properties param
%
%   F = UNCERTAIN(W,'name mixture', {param1}, {param2},...,mixtureweight)
%   declares W to be stochastic uncertainty with mixture distribution
%
%   INPUT
%    W : SDPVAR object
%
%   OUTPUT
%    F : Constraint object
%
%   EXAMPLE
%
%    Robust worst-case optimization
%
%    sdpvar x w
%    F = [x + w <= 1], W = [-0.5 <= w <= 0.5];
%    optimize([F,W,uncertain(w)],-x)
%
%    To specify random uncertainties, you specify the distribution, and all
%    distribution parameters following the syntax in the RANDOM command in
%    the Statistics Toolbox
%
%    sdpvar x w
%    F = [probability(x + w <= 1) >= 0.95, uncertain(w, 'uniform',0,1)];
%    optimize(F,-x);	
%
%    You can specify a function handle which generates samples to be used
%    with the SAMPLE command. YALMIP will always send a trailing argument
%    with dimensions.  
%
%    F = [x + w <= 1, uncertain(w,@mysampler,myarguments1,...)];
%
%    The standard uniform case above would thus be recovered with
%
%    F = [x + w <= 1, uncertain(w,@random,'uniform',0,1)];
%
%   See also PROBABILITY, SAMPLE

if nargin == 1 || ((nargin == 2) && strcmpi(varargin{1},'deterministic'))   
    x.extra.distribution.type = 'deterministic';    
else    
    if isa(varargin{1},'function_handle')
        temp = {varargin{:},x.dim};
    else
        temp = {@random,varargin{:},x.dim};
    end
    x.extra.distribution.type = 'stochastic';
    x.extra.distribution.generator = temp{1};
    x.extra.distribution.parameters = {temp{2:end-1}};
    x.extra.distribution.mixture = [];
    x.extra.distribution.characteristicfunction = [];
    x.extra.distribution.cdf = [];
    x.extra.distribution.icdf = [];
    x.extra.distribution.pcdf = [];
    
    if isequal(x.extra.distribution.generator, @random)
        % Check for a mixture definition
        if findstr('mix',x.extra.distribution.parameters{1})
            % Yes, mixture defined
            for i = 2:length(x.extra.distribution.parameters)
                if ~iscell(x.extra.distribution.parameters{i})
                    if i == length(x.extra.distribution.parameters)
                        % We support weights in both cell and vector, but
                        % temporarily place them in a cell for the
                        % dimension check below to be simple
                        x.extra.distribution.parameters{i} = num2cell(x.extra.distribution.parameters{i});
                    else
                        error('Mixture parameters should be placed in cells (except weights which can be a vector)')
                    end
                end
            end
            % Check weights
            mixtureweights = x.extra.distribution.parameters{end};            
            mixtureweights = cell2mat(mixtureweights);
            if ~all(abs(sum(mixtureweights,2)-1)<1e-12)
                error('Mixture weights should sum up to 1.')
            end
            if size(mixtureweights,1) < length(x.extra.distribution.parameters{2}{1})
                mixtureweights = repmat(mixtureweights, length(x.extra.distribution.parameters{2}{1}),1);
            end
            
            nMix = cellfun(@(c)size(c,2),x.extra.distribution.parameters);
            if ~all(nMix(2) == nMix(2:end))
                error('All parameter cells and mixture weights should have same length (#mixtures)')
            end
            
            % Remove mixture parameters and place in object instead
            x.extra.distribution.mixture = mixtureweights;
            x.extra.distribution.parameters = {x.extra.distribution.parameters{1:end-1}};
            x.extra.distribution.parameters{1} = strrep(x.extra.distribution.parameters{1},'mixture','');
            x.extra.distribution.parameters{1} = strrep(x.extra.distribution.parameters{1},'mix','');
            x.extra.distribution.parameters{1} = lower(strtrim(x.extra.distribution.parameters{1}));
            
            if strcmp(x.extra.distribution.parameters{1},{'mvnrnd','mvnrndfactor','dro','data','moment','momentf'})
                error(['Mixtures of ' x.extra.distribution.parameters{1} ' distributions not supported (yet...)']);
            end
                                    
            if isa(varargin{1},'function_handle')
                % Hmm, do we support generic mixtures of this type
                error('Mixture of sample generators not supported');
            else
                temp = {@random,x.extra.distribution.parameters{1},x.dim};
            end
        end        
    end
    try
        if any(cellfun('isclass',temp,'sdpvar')) || (strcmp(func2str(temp{1}),'random') && (any(strcmp(temp{2},{'mvnrnd','mvnrndfactor','dro','data','moment','momentf','laplace','cauchy'}))))
            % Don't try to evaluate special case distributions, such as
            % distributions with decision variables, or aditional cases
            % 'normalm' (multivariate normal) or 'normalf' (factor covar)
        else
            if isempty(x.extra.distribution.mixture)
                % Check if std. dev has been given as a matrix for a
                % normal. If so, it must be diagonal, otherwise mvnrnd
                % should be used
                if strcmp(x.extra.distribution.parameters{1},'normal')
                    stdDev = x.extra.distribution.parameters{3};
                    if min(size(stdDev)) > 1
                        msg1 = sprintf('Matrix standard deviation detected on ''normal'' declaration. A vector is expected\n');
                        msg2 = sprintf('''mvnrnd'' (with covariance parameter) should be used for truly multivariate Gaussians.\n');
                        msg3 = sprintf('''normal'' (with standard deviation parameter) is used for (possibly vectorized) scalar Gaussians');
                        msg = [msg1 msg2 msg3];
                        error(msg)
                    end
                end
                % Sanity checks done, now try to use this distribution
                % model by calling the data generator
                temp = feval(temp{:});
            else
                for i = 1:length(x.extra.distribution.mixture)
                    temp = x.extra.distribution.parameters;
                    for j = 2:length(x.extra.distribution.parameters)
                        temp{j} = temp{j}{i};
                        if isa(temp{j},'sdpvar')
                            temp{j} = value(temp{j});
                            temp{j}(isnan(temp{j})) = 0;
                        end
                    end
                    y = random(temp{:});
                end
            end
        end
    catch
        if isequal(temp{2},'normalm') || isequal(temp{2},'normalf')
            error('normalm and normalf have been replaced by mvnrnd and mvrndnfactor')
        else
            disp(lasterr);
            error('Trial evaluation of attached sample generator failed. Did you really specify reasonable parameters?')
        end
    end
    yalmip('addDistribution',  x, x.extra.distribution);
    x = lmi(x);
end