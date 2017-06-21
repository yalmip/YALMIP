function varargout=find(varargin)
%FIND (overloaded)
%
% [ind,val] = find(x,k,{'first','last'})
%
% Note: If x is zero, the index will be length(x)+1, and the value will be 0
% Hence, a way to force some element to be non-zero, is to add the
% constraint ind <= length(x).
% This also means that you never should use x(ind) if x can be all zeros,
% but use the second output val instread.
%
% FIND is implemented in the nonlinear operator framework using a big-M
% model.

switch class(varargin{1})
    
    case 'double'
        error('Overloaded SDPVAR/FIND CALLED WITH DOUBLE. Report error')
        
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        
        x = varargin{1};
        if min(size(x)) > 1
            error('SDPVAR/FIND currenty only support vector arguments');
        end
        x = reshape(x,1,[]);
         
        if nargin > 1
            if ~(isequal(varargin{2},[]) || isequal(varargin{2},1))
                error('SDPVAR/FIND currenty only support 1 non-zero element (second argument must be 1 or [])');
            end
        end
        k = 1;
        
        if nargin > 2
            if isequal(varargin{3},'last')
                error('SDPVAR/FIND currenty only supports detecting first non-zero element');
            end
        end
        which = 'first';
                            
        x_extended = [x 1];
        ind = yalmip('define',mfilename,x_extended,k,which);
        varargout{1} = ind;
        
        if nargout == 2
            x_0 = [x 0];
            varargout{2} = subsref(x_0,struct('type','()','subs',{{ind}}));
        end
        
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
