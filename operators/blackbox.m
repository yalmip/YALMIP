function varargout = blackbox(varargin)
%BLACKBOX Define user-specified function
%
% y = blackbox(x,f,properties)
%
% f     : Function handle 
% x     : SDPVAR
%
% Operator properties can be supplied to help the global 
% solver BMIBNB when brancing and performing bound propagation
%
% 'convex','concave','positive','negative','increasing','decreasing'
%
% Operator functions can also be supplied with the most important one
% being 'derivative' which is used also in local solvers.
%
% Example
%
% y = blackbox(x,@(x)sin(x).^2,'positive','derivative',@(x)2.*sin(x).*cos(x))
%
% See also SDPVAR

if ~isa(varargin{1},'char') && any(strcmp(cellfun(@class,varargin,'UniformOutput',0),'sdpvar'))
        
    dataExample =  zeros(size(varargin{1},1),size(varargin{1},2));%varargin_copy{1};
    try
        output = feval(varargin{2},dataExample);
    catch
        dataExample = value(varargin{1});    
        try
            output = feval(varargin{2},dataExample);
        catch
            warning('Cannot figure out output dimension. Assuming 1x1. Assign a reasonable value to the decision variable to avoid this warning');
            output = 0;
        end
    end
    
    % And now see what we got
    operator = CreateBasicOperator('callback');
    if nargin > 2
        reading = 3;
        while reading <= nargin
            switch class(varargin{reading})
                case 'function_handle'
                    error('Did not expect function handle here');
                case 'char'
                    switch varargin{reading}
                        case 'derivative'
                            operator.derivative = varargin{reading+1};
                            reading = reading + 2;
                        case {'convex','concave'}
                            operator.convexity = varargin{reading};
                            reading = reading + 1;
                        case {'positive','negative'}
                            operator.definiteness = varargin{reading};
                            reading = reading + 1;
                        case {'increasing','decreasing'}
                            operator.monotinicity = varargin{reading};
                            reading = reading + 1;
                        case 'bounds'
                            operator.bounds  = varargin{reading+1};
                            reading = reading + 2;
                        case 'convexhull'
                            operator.bound  = varargin{reading+1};
                            reading = reading + 2;
                        case 'stationary'
                            operator.stationary  = varargin{reading+1};
                            reading = reading + 2;                            
                        otherwise
                            error('Property not recognized')
                    end
                otherwise
            end
        end
    end
    % Populate completely to save for later 
    operator = assertOperatorProperties(operator);
    if numel(output) == 1
        if numel(dataExample)>1
            f = varargin{2};
            [n,m] = size(dataExample);
            g = @(x)f(reshape(x,n,m));
            varargin{2} = g;
        end        
        varargin = {varargin{1:2},operator};
        varargout{1} = yalmip('define',mfilename,varargin{:});
    else
        % MANY -> MANY
        y = [];
        [n,m] = size(varargin{1});
        varargin{1} = varargin{1}(:);
        for i = 1:length(varargin{1})
            varargin_copy = varargin;
            varargin_copy{1} = varargin_copy{1}(i);
            varargin_copy = {varargin_copy{1:2},operator};
            y = [y;yalmip('define',mfilename,varargin_copy{:})];
        end
        varargout{1} = reshape(y,n,m);
    end
    return  
end

switch class(varargin{1})
    
    case {'struct','double'} % What is the numerical value of this argument (needed for displays etc)
        varargout{1} = feval(varargin{2},varargin{1});
        
    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        
        operator = varargin{end};        
        if isempty(operator.bounds)
            canDeriveBoundsFromMonotone = strcmpi(operator.monotonicity,'increasing') || strcmpi(operator.monotonicity,'decreasing');            
            canDeriveBoundsFromConvex = strcmpi(operator.convexity,'convex') || strcmpi(operator.convexity,'concave');
            canDeriveBoundsFromConvex = canDeriveBoundsFromConvex && ~isempty(operator.stationary);
            if ~(canDeriveBoundsFromMonotone || canDeriveBoundsFromConvex)           
                operator.bounds = @(xL,xU)(bounds(xL,xU,varargin{4},operator));
                %operator.bounds = @(xL,xU)(bounds(xL,xU,varargin{4}));
            end
        end
        if isempty(operator.convexhull)
            canDeriveHullFromVexity = strcmpi(operator.convexity,'concave') || strcmpi(operator.convexity,'convex');
            canDeriveHullFromVexity = canDeriveHullFromVexity && ~isempty(operator.derivative);
            if ~canDeriveHullFromVexity
                operator.convexhull = @(xL,xU)(convexhull(xL,xU,varargin{4},operator));
                %operator.convexhull = @(xL,xU)(convexhull(xL,xU,varargin{4}));
            end
        end
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};
        
    otherwise
        error([upper(mfilename) ' called with weird argument']);
end

function [Ax,Ay,b,K] = convexhull(xL,xU,varargin)
% Note, properties are appended when blackbox
% operator is used. Should be exploited
if length(xL)==1 && ~(isinf(xL) || isinf(xU))
    N = 100;
    z = linspace(xL,xU,N);
    
    xS = [];   
    operator = varargin{end};
    if ~isempty(operator.stationary)
        % N.B stationarity does not mean min/max as 
        % we have no knowledge of convexity. If we had
        % this would have been exploited further up
        if isa(operator.stationary,'function_handle')
            % FIXME Add support
        else
            xS_ = operator.stationary;            
            loc = find((xS_ > xL) & (xS_ < xU));
            xS = xS_(loc);  
        end
    end
    z = sort([z xS]);
    
    fz = feval(varargin{1},z);
    % create 4 bounding planes
    % f(z) < k1*(x-XL) + f(xL)
    % f(z) > k2*(x-XL) + f(xL)
    % f(z) < k3*(x-XU) + f(xU)
    % f(z) > k4*(x-XU) + f(xU)
    % f(z) < k5*(x-XM) + f(xM)
    % f(z) > k6*(x-XM) + f(xM)
    
    k1 = max((fz(2:end)-fz(1))./(z(2:end)-xL))+1e-12;
    k2 = min((fz(2:end)-fz(1))./(z(2:end)-xL))-1e-12;
    k3 = min((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))+1e-12;
    k4 = max((fz(1:end-1)-fz(end))./(z(1:end-1)-xU))-1e-12;
    Ax = [-k1;k2;-k3;k4];
    Ay = [1;-1;1;-1];
    b =  [k1*(-z(1)) + fz(1);-(k2*(-z(1)) + fz(1));k3*(-z(end)) + fz(end);-(k4*(-z(end)) + fz(end))];
    K.f = 0;
    K.l = length(b);
else
    Ax = [];
    Ay = [];
    b = [];
    K = [];
end

function LU = bounds(xL,xU,varargin)
% Note, properties are appended when blackbox
% operator is used. Should be exploited
if length(xL)==1 && ~isinf(xL) & ~isinf(xU)
    operator = varargin{end};
    % FIXME: Add bisection, faster than sampling
    %if strcmp(operator.convexity,'convex')    
    %    LU = bisection_convex(xL,xU,varargin{1});
    %elseif strcmp(operator.convexity,'concave')
    %    LU = bisection_concave(xL,xU,varargin{1});
    %end
    N = 100;
    xtest = linspace(xL,xU,N);    
    xS = [];   
    if ~isempty(operator.stationary)
        % Stationarity does not mean min/max as 
        % we have no knowledge of convexity. If we had
        % this it would have been exploited further up
        if isa(operator.stationary,'function_handle')
            % FIXME Add support
        else
            xS_ = operator.stationary;
            loc = ((xS_ > xL) & (xS_ < xU));
            xS = xS_(loc);            
        end
    end    
    if xL<0 && xU > 0
        % Origin is often a critical point
        xS = [xS 0];        
    end
    xtest = [xtest xS];
    xtest = sort(xtest);
    values = feval(varargin{1},xtest);       
    [minval,minpos] = min(values);
    [maxval,maxpos] = max(values);
    % Concentrate search around the extreme regions
    xtestmin = linspace(xtest(max([1 minpos-5])),xtest(min([N minpos+5])),N);
    xtestmax = linspace(xtest(max([1 maxpos-5])),xtest(min([N maxpos+5])),N);
    xtestmin = [xtestmin xS];
    xtestmax = [xtestmax xS];
    values1 = feval(varargin{1},xtestmin);
    values2 = feval(varargin{1},xtestmax);
    L = min([values1 values2]);
    U = max([values1 values2]);
else
    L = -inf;
    U = inf;
end
LU = [L U];


function LU = bisection_convex(xL,xU,f)
fL = feval(f,xL);
fU = feval(f,xL);
LU(2) = max(fL,fU);

function LU = bisection_concave(xL,xU,f)
fL = feval(f,xL);
fU = feval(f,xL);
LU(1) = min(fL,fU);