function varargout=graph_square(varargin)
%GRAPH_SQUARE (overloaded)

% Used internally when users write stuff like norm(x,1)^2
switch class(varargin{1})
    case 'double'
        varargout{1} = varargin{1}.^2;
        
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.        
        varargout{1} = yalmip('define',mfilename,varargin{1});           

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
        switch varargin{1}
            case 'graph'
                % Description using epigraphs
                t = varargin{2};
                X = varargin{3};             
                varargout{1} = cone([t;X]);                
                varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','graph');
                varargout{3} = X;           
            otherwise
                error('SDPVAR/GRAPH_SQUARE called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/GRAPH_SQUARE');
end