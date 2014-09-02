function varargout = log(varargin)
%LOG (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/LOG CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
         % Try to detect logsumexp construction etc
        varargout{1} = check_for_special_cases(varargin{:});
        % Nope, then just define this logarithm
        if isempty(varargout{1})
            varargout{1} = InstantiateElementWise(mfilename,varargin{:});
        end
        
    case 'char'

        X = varargin{3};      
        F = (X >= 1e-8);

        operator = struct('convexity','concave','monotonicity','increasing','definiteness','none','model','callback');       
        operator.convexhull = @convexhull;
        operator.bounds = @bounds;
        operator.domain = [0 inf];
        operator.derivative = @(x)(1./(abs(x)+eps));
        operator.inverse = @(x)(exp(x));

        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/LOG called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
if xL <= 0
    % The variable is not bounded enough yet
    L = -inf;
else
    L = log(xL);
end
if xU < 0
    % This is an infeasible problem
    L = inf;
    U = -inf;
else
    U = log(xU);
end

function [Ax, Ay, b, K] = convexhull(xL,xU)
K = [];
if xL <= 0
    fL = inf;
else
    fL = log(xL);
end
fU = log(xU);
dfL = 1/(xL);
dfU = 1/(xU);
%xM = (xU - xL)/(fU-fL);
xM = (xL + xU)/2;
fM = log(xM);
dfM = 1/xM;

[Ax,Ay,b] = convexhullConcave(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
remove = isinf(b) | isinf(Ax) | isnan(b);
if any(remove)
    remove = find(remove);
    Ax(remove)=[];
    b(remove)=[];
    Ay(remove)=[];
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

