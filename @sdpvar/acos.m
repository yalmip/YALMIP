function varargout = acos(varargin)
%ACOS (overloaded)

% Author Johan Löfberg
% $Id: acos.m,v 1.6 2008-12-11 12:02:41 joloef Exp $
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ACOS CALLED WITH DOUBLE. Report error')

    case 'sdpvar'
        % We will compute value in an internal file in order to extend
        % the function beyond -1 and 1 (since fmincon makes calls outside
        % bounds...)
        varargout{1} = InstantiateElementWise('acos_internal',varargin{:});

%     case 'char'
% 
%         operator = struct('convexity','none','monotonicity','decreasing','definiteness','none','model','callback');
%         operator.convexhull = [];
%         operator.bounds = @bounds;
%         operator.derivative = @derivative;
%         operator.range = [-pi/2 pi/2];
%         operator.domain = [-1 1];
% 
%         varargout{1} = [];
%         varargout{2} = operator;
%         varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/ACOS called with CHAR argument?');
end
% 
% function [L,U] = bounds(xL,xU)
% L = real(acos(xU));
% U = real(acos(xL));
% 
% function df = derivative(x)
% 
% df = (-(1 - x.^2).^-0.5);
% % FMINCON is really stupid sometimes and make calls outside domain!
% % df(x<-1) = -1e4;
% % df(x>1) = -1e4;
