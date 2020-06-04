function varargout = eig_yalmip_internal(varargin)

switch class(varargin{1})

    case 'double'
        X = varargin{1}(:);
        n = sqrt(length(X));
        X = reshape(X,n,n);
        varargout{1} = eig((X+X')/2);

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.

        x = varargin{1};
        y = yalmip('definemulti',mfilename,x,[size(x,1) 1]);
        varargout{1} = y;

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
        %   disp('The EIG operator is not supported in optimization problems.')
        %   error

        operator = CreateBasicOperator('callback');
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

        % Inofficial way to model several nonlinear variables in
        % one call
        % varargout{2}.models = getvariables(t):getvariables(t)+length(X)-1;

    otherwise
        error([upper(mfilename) ' called with weird argument']);
end
