function varargout=quadratic_over_affine_expanded(varargin)
%quadratic_over_affine (internal)

switch class(varargin{1})
    case 'double'
        varargout{1} = varargin{1}'*varargin{1}/varargin{2};
        
    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
        switch varargin{1}
            case 'graph'
                % Description using epigraphs
                t = varargin{2};
                q = varargin{3};
                y = varargin{4};
                varargout{1} = [cone([2*q;t-y],t+y)];
                varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','graph');
                varargout{3} = [q;y];
              
            case {'exact','integer','callback'}
                
                t = varargin{2};
                q = varargin{3};
                y = varargin{4};
                varargout{1} = [];
                varargout{2} = struct('convexity','none','monotonicity','none','definiteness','positive','model','callback');
                varargout{3} = [q;y];

            otherwise
                error('SDPVAR/ABS called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/ABS');
end
