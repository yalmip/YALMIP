function varargout = SIGN(varargin)
%SIGN (overloaded)

% Author Johan Löfberg
% $Id: sign.m,v 1.4 2009-02-12 10:31:47 joloef Exp $
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/SIGN CALLED WITH DOUBLE. Report error')

    case 'sdpvar' 
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char' 
        switch varargin{1}
            case {'graph','exact'}

                t = varargin{2};
                X = varargin{3};

                d = binvar(1,1);
                [M,m] = derivebounds(X);
                F = [X >= d*m,X <=(1-d)*M, t == 1-2*d];

                varargout{1} = F;
                varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none');
                varargout{3} = X;

            otherwise
                error('SDPVAR/SIGN called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/SIGN');
end
