function varargout = invpos(varargin)
% INVPOS  Returns model of 1./x

% Author Johan Löfberg
% $Id: invpos.m,v 1.1 2007-12-17 12:16:41 joloef Exp $

switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        varargout{1} = 1./varargin{1};

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});
        
    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.
            
            varargout{1} = set(cone([2;X-t],X+t));
            varargout{2} = struct('convexity','convex','monotonicity','decreasing','definiteness','positive','model','graph');
            varargout{3} = X;
        else
        end
    otherwise
end
