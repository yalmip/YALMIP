function varargout = fraction_internal2(varargin)

% Author Johan Löfberg
% $Id: fraction_internal2.m,v 1.1 2007-08-22 11:31:19 joloef Exp $
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
