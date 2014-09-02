function varargout = gammainc_x(varargin)
%GAMMAINC_X

switch class(varargin{1})

    case 'double'       
        varargout{1} = gammainc(varargin{1},varargin{2});
  
    case 'char'

        operator = struct('convexity','none','monotonicity','none','definiteness','none','model','callback');
        operator.convexhull = [];              
        operator.range = [0 1];
        operator.domain = [1e-6 inf];
        
        varargout{1} = [];
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('SDPVAR/GAMMAINC_X called with strange argument?');
end
