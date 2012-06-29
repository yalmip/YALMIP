function varargout = openopt_fun(varargin)

persistent params

if nargin>1
    params = varargin{2};
    return
end

if ischar(varargin{1})
    switch varargin{1}
        case 'getStartPoint'
            varargout{1} = params.x0;
            return
        case 'getOptimPoint'
            varargout{1} =[];
            return
    end
end

xevaled = zeros(1,length(params.c));
xevaled(params.linearindicies) = varargin{1};

% Apply the precomputed evaluation scheme (if necessary)
xevaled = apply_recursive_evaluation(params,xevaled);

xevaled = xevaled(:);
varargout{1} = params.f + (params.c'+xevaled'*params.Q)*xevaled;
