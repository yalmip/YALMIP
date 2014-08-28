function varargout = optimizer_operator(varargin)

switch class(varargin{1})

    case 'double'
        x = varargin{1}(:);
        O = varargin{2};
        [varargout{1},problem] = O{x};
        if problem
            varargout{1} = inf(struct(O).dimout);
        end

    case {'sdpvar','optimizer'} % Overloaded operator for SDPVAR objects.
        O = struct(varargin{1});
        varargout{1} = yalmip('definemulti',mfilename,varargin{2},varargin{1},O.dimin);

    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph

        X = varargin{3};

        % Let YALMIP know about convexity etc
        varargout{1} = [];
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        varargout{3} = X;

    otherwise
        error('optimizer_operator called with CHAR argument?');
end