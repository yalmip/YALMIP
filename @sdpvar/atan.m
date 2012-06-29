function varargout = atan(varargin)
%ATAN (overloaded)

% Author Johan Löfberg
% $Id: atan.m,v 1.8 2007-08-02 18:16:26 joloef Exp $
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ATAN CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        operator = struct('convexity','none','monotonicity','increasing','definiteness','none','model','callback');
        operator.convexhull = [];
        operator.bounds = @bounds;
        operator.derivative = @(x)((1+x.^2).^-1);
        operator.range = [-pi/2 pi/2];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ATAN called with CHAR argument?');
end

function [L,U] = bounds(xL,xU)
L = atan(xL);
U = atan(xU);