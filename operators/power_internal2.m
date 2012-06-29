function varargout = power_internal2(varargin)
%power_internal2 (overloaded)
% Used for cases such as x^x and 2^x, which is treated as evaluation-based
% operators

% Author Johan Löfberg
% $Id: power_internal2.m,v 1.7 2007-08-02 19:17:36 joloef Exp $
switch class(varargin{1})

    case 'double'
        varargout{1} = real(varargin{1}.^varargin{2});

    case 'sdpvar' %
        if isa(varargin{2},'sdpvar')
            error('x^y currently not supported for SDPVAR x and SDPVAR y')
        else
            varargout{1} = yalmip('define',mfilename,varargin{:});
        end

    case 'char'

        X = varargin{3};
        Y = varargin{4};

        varargout{1} = [];
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','exact');
        varargout{3} = [X(:);Y(:)];
    otherwise
        error('SDPVAR/power_interna2 called with CHAR argument?');
end
