function varargout = interp1(varargin)
%INTERP1 (overloaded)

switch class(varargin{3})

    case 'sdpvar'   
        
        if numel(varargin{3})>1
            error('Only scalar functions supported in SDPVAR/INTERP1');
        end
        
        if nargin == 3
            varargin{4} = 'linear';
        end
        
        if ~isa(varargin{1},'double') || ~isa(varargin{2},'double')
            error('First 2 arguments in interp1 should be approximation data');
        end
        
        if any(diff(varargin{1})<0)
            error('First arguments in interp1 should be monotonically increasing');
        end
        
        if isequal(varargin{4},'graph')
            if  ~isconvexdata(varargin{1},varargin{2}) && ~isconvexdata(varargin{1},-varargin{2})
                error('Data has to be convex or concave for graph approximant');
            end
        end
        
        % Reorder arguments to make sure sdpvar is first argument.
        % Use local version of interp to deal with this
        % Arguments xi yi x method -> x xi yi method
        varargout{1} = yalmip('define','interp1_internal',varargin{[3 1 2 4]});   
     
    otherwise
        error('SDPVAR/INTERP1 called with strange argument!');
end

function isconvex = isconvexdata(xi,yi)
finterp = (yi(1:end-2) + yi(3:end))/2;
if all(finterp >= yi(2:end-1))
    isconvex = 1;
else
    isconvex = 0;
end