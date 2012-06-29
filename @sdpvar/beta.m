function varargout=beta(varargin)
%BETA (overloaded)

% Author Johan Löfberg
% $Id: beta.m,v 1.4 2008-05-01 21:57:37 joloef Exp $
switch class(varargin{1})

    case 'sdpvar'
        if ~isa(varargin{2},'double')
            error('W is not allowed to be an SDPVAR object')
        end
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char'

        X = varargin{3};
        F = set(X >= eps);
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