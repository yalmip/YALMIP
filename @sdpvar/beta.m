function varargout=beta(varargin)
%BETA (overloaded)

switch class(varargin{1})

    case 'sdpvar'
        if ~isa(varargin{2},'double')
            error('W is not allowed to be an SDPVAR object')
        end
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        X = varargin{3};
        F = (X >= eps);
        operator = struct('convexity','convex',...
                          'monotonicity','decreasing',...
                          'definiteness','positive',...
                          'model','callback');
        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = varargin{3};

    otherwise
        error('Strange type on first argument in SDPVAR/BETA');
end