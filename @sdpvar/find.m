function varargout=find(varargin)
%FIND (overloaded)
%
% ind = sort(x,k,{'first','last'})
%
% FIND is implemented in the nonlinear operator framework using a big-M
% model.

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/FIND CALLED WITH DOUBLE. Report error')

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        
        if min(size(varargin{1})) > 1
            error('SDPVAR/FIND currenty only support vector arguments');
        end
        
        if nargin > 1
            if ~(isequal(varargin{2},[]) || isequal(varargin{2},1))
                error('SDPVAR/FIND currenty only support 1 non-zero element (second argument must be 1 or [])');
            end
        end
                
        if nargin > 2
            if isequal(varargin{3},'last')                
                error('SDPVAR/FIND currenty only supports detecting first non-zero element');
            end
        end
        
        varargin{1} = reshape(varargin{1},1,[]);
        
          
        x = varargin{1};                   
        ind = yalmip('define',mfilename,x);
        varargout{1} = ind;        

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph

        t = varargin{2};
        X = varargin{3};
       
        % Call external to allow subsrefs in classs        
        F = find_internal(t,X);

        varargout{1} = F;
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = X;
    otherwise
        error('Strange type on first argument in SDPVAR/SORT');
end
