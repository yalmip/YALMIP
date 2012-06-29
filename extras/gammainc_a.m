function varargout = gammainc_a(varargin)
%GAMMAINC_A

% Author Johan Löfberg
% $Id: gammainc_a.m,v 1.8 2007-08-02 18:16:26 joloef Exp $
switch class(varargin{1})

    case 'double'       
        varargout{1} = gammainc(varargin{2},varargin{1});
  
    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.convexhull = [];     
        operator.range = [0 1];
        operator.domain = [1e-6 inf];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/GAMMAINC_A called with strange argument?');
end
