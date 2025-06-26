function varargout=std(varargin)
%STD (overloaded)
%
% V = std(x)

switch class(varargin{1})

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        X = varargin{1};
        if nargin > 1 || min(size(X))>1
            error('SDPVAR/STD only supports simple 1-D variance'),
        end
        varargout{1} = yalmip('define',mfilename,varargin{1});
    case 'char'
        switch varargin{1}
            case 'graph'
                t = varargin{2};
                X = varargin{3};
                F = cone([t;X])
                varargout{1} = F;
                varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','graph');
                varargout{3} = X;

            case 'exact'
                t = varargin{2};
                X = varargin{3};
                e = sdpvar(length(X),1);
                F = [e'*e == t^2, e == X, t>=0];
                varargout{1} = F;
                varargout{2} = struct('convexity','convex','definiteness','positive','model','exact');
                varargout{3} = X;
            otherwise
                error('SDPVAR/NORM called with weird argument?');
        end
    otherwise
        error('Weird first argument in SDPVAR/STD')
end
