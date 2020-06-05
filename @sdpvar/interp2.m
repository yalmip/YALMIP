function varargout = interp2(varargin)

switch class(varargin{4})
    
    case 'sdpvar'
        if nargin == 5
            varargin{6} = 'linear';
        end
        
        if isequal(varargin{6},'graph')
            if  ~isconvexmeshdata(varargin{1},varargin{2},varargin{3})
                error('Triangulation of data has to be convex or concave for graph approximant');
            end
        end
        
        if any(any(isnan(varargin{1}))) ||  any(any(isnan(varargin{2})))
            error('Interpolation grid contains NaNs');
        end
        if any(any(isnan(varargin{3})))
            error('Interpolation data contains NaNs');
        end
        
        
        % Reorder arguments to make sure sdpvar is first argument.
        % Use internal version of interp2 to deal with this
        % Arguments xi yi zi xv yv method -> [xv;yv] xi yi zi method
        decisionVariable = [varargin{4};varargin{5}];
        reordered = {decisionVariable,varargin{1:3},varargin{6}};
        varargout{1} = yalmip('define','interp2_internal',reordered{:});
        
    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end