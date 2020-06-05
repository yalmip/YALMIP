function varargout = eig(varargin)

switch class(varargin{1})

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.

        x = varargin{1};
        y = yalmip('definemulti','eig_yalmip_internal',x,[size(x,1) 1]);
        varargout{1} = y;

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
                     
        operator = CreateBasicOperator('callback');

        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

        % Inofficial way to model several nonlinear variables in
        % one call
        varargout{2}.models = getvariables(t):getvariables(t)+length(X)-1;

    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end
