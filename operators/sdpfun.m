function varargout = sdpfun(varargin)
%SDPFUN Gateway to general (elementwise) functions on SDPVAR variables (overloaded)

if ~isa(varargin{1},'char') && any(strcmp(cellfun(@class,varargin,'UniformOutput',0),'sdpvar'))
   
    fun_handles = zeros(nargin,1);
    for i = 1:nargin
        if isstr(varargin{end}) && isequal(varargin{i}(1),'@')
            fun = eval(varargin{i});
            varargin{i} = fun;
            fun_handles(i) = 1;
        elseif isa(varargin{i},'function_handle')
            fun_handles(i) = 1;
        end
    end
    % Temporarily remove derivative info (to avoid  change of legacy
    % code below)
    if nnz(fun_handles) == 2
        derivative = varargin{end};
        varargin = {varargin{1:end-1}}
    else
        derivative = [];
    end
    for i = 1:length(varargin)
        % we don't know in which order arguments are applied. If some argument
        % is an sdpvar, it means the user is defining the operator
        if isa(varargin{i},'sdpvar')
            % Overloaded operator for SDPVAR objects. Pass on args and save them.
            % try to figure out size of expected output (many->1 or many->many          
            varargin_copy = varargin;            
            for j = 1:length(varargin)-1
                % Create zero arguments
                if isa(varargin_copy{j},'sdpvar')
                    varargin_copy{j} = zeros(size(varargin_copy{j},1),size(varargin_copy{j},2));
                    dataExample = varargin_copy{j};
                end
            end
            output = feval(varargin_copy{end},varargin_copy{1:end-1});
            if prod(size(output)) == 1
                % MANY -> 1
                % Append derivative (might be empty)
                varargin{end+1} = derivative;
                if numel(dataExample)>1
                    f = varargin{end-1};
                    [n,m] = size(dataExample);
                    g = @(x)f(reshape(x,n,m));
                    varargin{end-1} = g;
                end
                varargout{1} = yalmip('define',mfilename,varargin{:});
            else
                % MANY -> MANY
                y = [];
                [n,m] = size(varargin{1});
                varargin{1} = varargin{1}(:);
                for i = 1:length(varargin{1})
                    varargin_copy = varargin;
                    varargin_copy{1} = varargin_copy{1}(i);
                    varargin_copy{end+1} = derivative;
                    y = [y;yalmip('define',mfilename,varargin_copy{:})];
                end
                varargout{1} = reshape(y,n,m);
            end
            return
        end
    end
end

switch class(varargin{1})

    case {'struct','double'} % What is the numerical value of this argument (needed for displays etc)
        varargout{1} = feval(varargin{end-1},varargin{1:end-2});

    case 'sdpvar' 
       % Should not happen
       error('Report bug. Call to SDPFUN with SDPVAR object');
       
    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
      
        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.bounds = @bounds;
        operator.convexhull = @convexhull;
        if ~isempty(varargin{end})
            operator.derivative = varargin{end};
        end

        varargout{1} = [];
        varargout{2} = operator;
        % Figure out which argument actually is the SDPVAR object (not
        % necessarily the first argument
        for i = 3:length(varargin)
            if isa(varargin{i},'sdpvar')
                varargout{3} = varargin{i};
                return
            end
        end
        error('SDPFUN arguments seem weird. Could not find any SDPVAR object');
        

    otherwise
        error('SDPVAR/SDPFUN called with CHAR argument?');
end

function [Ax,Ay,b] = convexhull(xL,xU,varargin)

if length(xL)==1 && ~(isinf(xL) | isinf(xU))
    z = linspace(xL,xU,100);
    fz = feval(varargin{end-1},z,varargin{1:end-2});
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
else
    Ax = [];
    Ay = [];
    b = [];
end

function [L,U] = bounds(xL,xU,varargin)

if length(xL)==1 && ~isinf(xL) & ~isinf(xU) 
    xtest = linspace(xL,xU,100);
    values = feval(varargin{end-1},xtest,varargin{1:end-2});
    [minval,minpos] = min(values);
    [maxval,maxpos] = max(values);
    % Concetrate search around the extreme regions
    xtestmin = linspace(xtest(max([1 minpos-5])),xtest(min([100 minpos+5])),100);
    xtestmax = linspace(xtest(max([1 maxpos-5])),xtest(min([100 maxpos+5])),100);
    values1 = feval(varargin{end-1},xtestmin,varargin{1:end-2});
    values2 = feval(varargin{end-1},xtestmax,varargin{1:end-2});
    L = min([values1 values2]);
    U = max([values1 values2]);
else
    L = -inf;
    U = inf;
end
