function varargout = log(varargin)

switch class(varargin{1})

    case 'sdpvar'
         % Try to detect logsumexp construction etc
        varargout{1} = check_for_special_cases(varargin{:});
        % Nope, then just define this logarithm
        if isempty(varargout{1})
            varargout{1} = InstantiateElementWise(mfilename,varargin{:});
        end
        
    case 'char'
                    
        operator = CreateBasicOperator('concave','increasing','callback');                             
        operator.derivative = @(x)(1./(abs(x)+eps));
        operator.inverse = @(x)(exp(x));
        operator.domain = [0 inf];
        operator.singularity = [0 -inf -inf] ;

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end

function f = check_for_special_cases(x)
f = [];
% Check for log(1+x)
base = getbase(x);
if all(base(:,1)==1)
    f = slog(x-1);
    return;
end
% Check if user is constructing log(sum(exp(x)))
if base(1)~=0
    return
end
if ~all(base(2:end)==1)
    return
end
modelst = yalmip('extstruct',getvariables(x));
if isempty(modelst)
    return;
end
if length(modelst)==1
    models{1} = modelst;
else
    models = modelst;
end
% LOG(DET(X))
if length(models)==1
    if strcmp(models{1}.fcn,'det_internal')
        n = length(models{1}.arg{1});
        try
            f = logdet(reshape(models{1}.arg{1},sqrt(n),sqrt(n)));
        catch
        end
        return
    end
end
% LOG(EXP(x1)+...+EXP(xn))
for i = 1:length(models)
    if ~strcmp(models{i}.fcn,'exp')        
        return
    end      
end
p = [];
for i = 1:length(models)
    p = [p;models{i}.arg{1}];
end
f = logsumexp(p);

