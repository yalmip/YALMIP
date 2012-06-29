function varargout=ceil(varargin)
%CEIL (overloaded)

% Author Johan Löfberg
% $Id: ceil.m,v 1.4 2007-07-26 17:10:13 joloef Exp $

switch class(varargin{1})
    
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        
        x = varargin{1};
        dim = size(x);
        x = reshape(x,prod(dim),1);
        y = [];
        for i = 1:prod(dim)
            y = [y;yalmip('define',mfilename,extsubsref(x,i))];
        end
        y = reshape(y,dim);
        varargout{1} = y;
        
    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
        switch varargin{1}
            case {'milp','graph'}
                % Description using epigraphs
                t = varargin{2};
                X = varargin{3};
                
                F = set([X <= t <= X + 1]) + set(integer(t));
                
                varargout{1} = F;
                varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
                varargout{3} = X;
                
            otherwise
                error('SDPVAR/CEIL called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/CEIL');
end
