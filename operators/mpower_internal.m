function varargout = mpower_internal(varargin)
switch class(varargin{1})
    case 'double'
        varargout{1} = varargin{1};
    case 'char'
        % Note, the order here is assumed in GP formulation
        varargout{1} = (varargin{3} == varargin{2});        
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','exact');
        varargout{3} = varargin{3};
    otherwise
end
