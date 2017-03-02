function varargout = interp2(varargin)
%INTERP2 (overloaded)

switch class(varargin{4})

    case 'sdpvar'                         
        if nargin == 5
            varargin{6} = 'linear';
        end        
        % Reorder arguments to make sure sdpvar is first argument.
        % Use local version of interp to deal with this
        % Arguments xi yi zi xv yv method -> [xv;yv] xi yi zi method
        decisionVariable = [varargin{4};varargin{5}];
        reordered = {decisionVariable,varargin{1:3},varargin{6}};
        varargout{1} = yalmip('define','interp2_internal',reordered{:});   
     
    otherwise
        error('SDPVAR/INTERP called with strange argument!');
end