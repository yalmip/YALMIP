function x = uncertain(x,varargin)
%UNCERTAIN Declares a variable as uncertain
%
%   F = UNCERTAIN(W) is used to describe the set of uncertain variables
%   in an uncertain program, both deterministic and stochastic
%
%   INPUT
%    W : SDPVAR object
%    S : Optional distribution information for random uncertainty
%
%   OUTPUT
%    F : Constraint object
%
%   Uncertain is used to declare uncertain variables in robust
%   deterministic worst-case optimization. It can also be used to specify
%   uncertain random variables, and their associated distribution.
%
%   EXAMPLE
%
%    Robust worst-case optimization
%
%    sdpvar x w
%    F = [x + w <= 1], W = [-0.5 <= w <= 0.5];
%    optimize([F,W,uncertain(w)],-x)
%
%   To specify random uncertainties, you specify the distribution, and all
%   distribution parameters following the syntax n the RANDOM command in
%   the Statistics Toolbox
%
%    sdpvar x w
%    F = [x + w <= 1, uncertain(w, 'uniform',0,1)];
%    P = optimizer([F,W,uncertain(w)],-x,[],w,x)
%    S = sample(P,10); % Sample ten instances and concatenate models
%    S([])             % Solve and return optimal x
%
%    Alternatively, you can specify a function handle which generates
%    samples. YALMIP will always send a trailing argument with dimensions
%
%    F = [x + w <= 1, uncertain(w,@mysampler,myarguments1,...)];
%
%    The standard uniform case above would thus be recovered with
%
%    F = [x + w <= 1, uncertain(w,@random,'uniform',0,1)];
%
%
%   See also OPTIMIZE, ROBUSTMODEL, OPTIMIZER, SAMPLE

if nargin == 1 || ((nargin == 2) && strcmpi(varargin{1},'deterministic'))
    %  x.typeflag = 15;
    x.extra.distribution.type = 'deterministic';
    %  x = lmi(x);
else
    %  x.typeflag = 16;
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
    
    if isequal(x.extra.distribution.generator ,@random)
        % Check for a mixture definition
        if findstr('mix',x.extra.distribution.parameters{1})
            % Yes, mixture defined
            for i = 2:length(x.extra.distribution.parameters)
                if ~iscell(x.extra.distribution.parameters{i})
                    error('Mixture parameters should be placed in cells, including trailing mixture weights.')
                end
            end
            % Remove mixture parameters and place in object instead
            x.extra.distribution.mixture = [x.extra.distribution.parameters{end}{:}];
            x.extra.distribution.parameters = {x.extra.distribution.parameters{1:end-1}};
            x.extra.distribution.parameters{1} = strrep(x.extra.distribution.parameters{1},'mixture','');
            x.extra.distribution.parameters{1} = strrep(x.extra.distribution.parameters{1},'mix','');
            x.extra.distribution.parameters{1} = lower(strtrim(x.extra.distribution.parameters{1}));
            
            if isa(varargin{1},'function_handle')
                % Hmm, do we support generic mixtures of this type
                error
            else
                temp = {@random,x.extra.distribution.parameters{1},x.dim};
            end
        end
        switch x.extra.distribution.parameters{1}
            case 'exponential'
                % N.B parameterized in scale parameter beta = 1/lambda in
                % Statistics toolbox so we use that here
                x.extra.distribution.characteristicfunction = @(t,beta)1./(1-t*1i.*beta(:));
                x.extra.distribution.characteristicfunction_derivative = @(t,beta)1i*(1./beta(:))./(1./beta(:)-1i*t).^2;
            case 'laplace'
                x.extra.distribution.characteristicfunction = @(t,mu,b)exp(1i*mu(:)*t)./(1+b(:).^2.*t.^2);
                x.extra.distribution.characteristicfunction_derivative = [];
            case 'logistic'
                x.extra.distribution.characteristicfunction = @(t,mu,s)(exp(1i*mu(:).*t))./guarded_sinhc(pi*s(:).*t);
                x.extra.distribution.characteristicfunction_derivative = @(t,mu,s)(1i*mu(:).*phi(t) + exp(1i*mu(:)*t).*(pi*s(:).*sinh(pi*s(:)*t)-pi^2*s(:).^2.*cosh(pi*s(:)*t))./sinh(pi*s(:)*t).^2);
            case 'uniform'
                x.extra.distribution.characteristicfunction = @(t,a,b)(guarded_expdiv(b(:).*t,a(:).*t,t.*(b-a)));
                x.extra.distribution.characteristicfunction_derivative = @(t,a,b) (1 ./ (1i*(b(:)-a(:)) .* t.^2)) .* (t*(1i*b(:).*exp(1i * b(:) .* t) - 1i*a(:) .* exp(1i * a(:) .* t)) - (exp(1i * b(:) * t) - exp(1i * a(:) * t)));
                
            otherwise
                x.extra.distribution.characteristicfunction = [];
                x.extra.distribution.characteristicfunction_derivative = [];
        end
    end
    try
        if any(cellfun('isclass',temp,'sdpvar')) || (strcmp(func2str(temp{1}),'random') && (any(strcmp(temp{2},{'dro','data','moment','momentf','normalm','normalf','laplace'}))))
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
        disp(lasterr);
        error('Trial evaluation of attached sample generator failed. Did you really specify correct parameters?')
    end
    yalmip('addDistribution',  x, x.extra.distribution);
    %  x = lmi(x);
end

function y = guarded_sinhc(x)
y = sinh(x)./x;y(x==0)=1;

function y = guarded_expdiv(z1,z2,z3)
y = (exp(1i*z1)-exp(1i*z2))./(1i*z3);
y(z3==0) = 1;
