function varargout = mvncdf(varargin)
%MVNCDF
%
% y = MVNCDF(x)
%
% Computes/declares multivariate normal cumulative distribution function

switch class(varargin{1})

    case {'sdpvar','ndsdpvar'}
 
        if nargin > 1
            error('sdpvar/mvncdf currently only supports 1 argument, i.e. assumed zero mean and unit variance');
        end
        varargin{1} = reshape(varargin{1},[],1);
        varargout{1} = yalmip('define',mfilename,varargin{1});        

    case 'char'

        X = varargin{3};
        F = (X >= 0);

        operator = struct('convexity','none','monotonicity','none','definiteness','positive','model','callback');
        operator.range = [0 1];
        operator.domain = [0 inf];                        
        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/MVNCDF called with CHAR argument?');
end