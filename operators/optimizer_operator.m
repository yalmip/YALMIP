function varargout = optimizer_operator(varargin)

% Author Johan Löfberg
% $Id: optimizer_operator.m,v 1.6 2007-08-03 11:17:33 joloef Exp $

switch class(varargin{1})

    case 'double'
        x = varargin{1}(:);
        O = varargin{2};
        [varargout{1},problem] = O{x};
        if problem
            varargout{1} = inf;
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