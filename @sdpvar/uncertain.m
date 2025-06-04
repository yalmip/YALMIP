function x = uncertain(x,varargin)
%UNCERTAIN Declares a variable as uncertain
%
%   F = UNCERTAIN(W) declares W to be (deterministic) uncertainty
%
%   F = UNCERTAIN(W,name, param1, param2,...) declares W to be
%   stochastic uncertainty with distribution properties
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
                    error('Mixture parameters should be placed in cells, including trailing mixture weights.')
                end
            end
            % Check weights
            alpha = varargin{end};            
            if ~(abs(sum([alpha{:}])-1)<1e-12)
                error('Mixture weights should sum up to 1.')
            end
            
            nMix = cellfun(@length,x.extra.distribution.parameters);
            if ~all(nMix(2) == nMix(2:end))
                error('All parameter cells in mixture should have same length (#mixtures)')
            end
            
            % Remove mixture parameters and place in object instead
            x.extra.distribution.mixture = [x.extra.distribution.parameters{end}{:}];
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
                temp = feval(temp{:});
            else
                for i = 1:length(x.extra.distribution.mixture)
                    temp = x.extra.distribution.parameters;
                    for j = 2:length(x.extra.distribution.parameters)
                        temp{j} = temp{j}{i};
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
            error('Trial evaluation of attached sample generator failed. Did you really specify correct parameters?')
        end
    end
    yalmip('addDistribution',  x, x.extra.distribution);
    %  x = lmi(x);
end