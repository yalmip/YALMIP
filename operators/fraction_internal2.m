function varargout = fraction_internal2(varargin)

switch class(varargin{1})

    case 'double'
        varargout{1} = varargin{1}(1)./varargin{1}(2);

    case 'sdpvar' %
        if isa(varargin{2},'sdpvar')
            error('x^y currently not supported for SDPVAR x and SDPVAR y')
        else
            varargout{1} = yalmip('define',mfilename,varargin{:});
        end

    case 'char'

        X = varargin{3};

        varargout{1} = [];
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','exact');
        varargout{3} = [X(:)];
    otherwise
        error('SDPVAR/fraction_internal2 called with CHAR argument?');
end
