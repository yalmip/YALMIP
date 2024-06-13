function varargout=rem(varargin)
%REM (overloaded)

switch class(varargin{1})

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.

        if ~isa(varargin{2},'double')
            error('REM is currently only supported for DOUBLE second argument');
        end

        x = varargin{1};
        y = varargin{2};
        % Some boring code for scalarization
        if prod(size(x))==1 & prod(size(y))>1
            x = x*ones(size(y));
        elseif prod(size(y))==1 & prod(size(x))>1
            y = y*ones(size(x));
        end
        if ~all(size(x) == size(y))
            error('Matrix dimensions must agree.');
        end
        dim = size(x);
        x = reshape(x,prod(dim),1);
        y = reshape(y,prod(dim),1);
        z = [];
        % Create one variable for each element
        for i = 1:length(x)
            xi = extsubsref(x,i);
            yi = extsubsref(y,i);
            inarg = {xi,yi};
            z = [z;yalmip('define',mfilename,inarg{:})];
        end
        z = reshape(z,dim);
        varargout{1} = z;

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph

        % Description using epigraphs
        t = varargin{2};
        x = varargin{3};
        y = varargin{4};

        % t = rem(x,y), i.e. t = x - n*y, n = fix(x/y)
        d1 = binvar(1,1);
        d2 = binvar(1,1);
        [M,m] = derivebounds(x);
        n = intvar(1,1);
        F = [t == x - y*n, d1+d2 == 1];
        F = [F,x>=m*(1-d1), x<=M*(1-d2),(x/y)-d1 <= n <= (x/y)+d2];
        if ~isinf(m) && isa(y,'double')
            F = [F, fix(m/y) <= n];
        end
        if ~isinf(M) && isa(y,'double')
            F = [F, n <= fix(M/y)];
        end

        varargout{1} = F;
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = [x(:);y(:)];

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end
