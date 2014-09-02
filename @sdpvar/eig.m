function varargout = eig(varargin)
%EIG (overloaded)

switch class(varargin{1})

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.

        x = varargin{1};
        y = yalmip('definemulti','eig_yalmip_internal',x,[size(x,1) 1]);
        varargout{1} = y;

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
     
        X = varargin{3};

        varargout{1} = [];
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        varargout{3} = X;

        % Inofficial way to model several nonlinear variables in
        % one call
        varargout{2}.models = getvariables(t):getvariables(t)+length(X)-1;

    otherwise
        error('Strange type on first argument in SDPVAR/EIG');
end
