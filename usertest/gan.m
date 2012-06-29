function varargout = gan(varargin)
% Author Johan Löfberg
% $Id: gan.m,v 1.6 2007/08/21 14:49:42 joloef Exp $
switch class(varargin{1})

    case 'double'
        x = varargin{1};
        a = varargin{2};
        varargout{1} = sum(x(:).*a(:).^(1./x(:)));

    case 'sdpvar'
        X = varargin{1};
        if min(size(X))>1
            error('ENTROPY only defined for vector arguments');
        else
            varargout{1} = yalmip('define',mfilename,varargin{:});
        end

    case 'char'

        X = varargin{3};
        F = set(X > 0);

        varargout{1} = F;
        varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','none','model','callback');
        varargout{3} = X;

    otherwise
        error('SDPVAR/LOG called with CHAR argument?');
end
