function varargout = pow10(varargin)
%POW10 (overloaded)

switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        % SHOULD NEVER HAPPEN, THIS SHOULD BE CAUGHT BY BUILT-IN
        varargout = 10.^varargin{1};

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        if length(varargin{1}) == 1
            varargout{1} = yalmip('addEvalVariable',mfilename,varargin{1});
        else
            y = [];
            for i = 1:length(varargin{1})
                y = [y;yalmip('addEvalVariable',mfilename,extsubsref(varargin{1},i))];
            end
            varargout{1} = y;
        end

    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        switch varargin{1}
            case 'graph'
                t = varargin{2};
                X = varargin{3};
                
                % This is different from so called extended operators
                % Just do it!
                F = SetupEvaluationVariable(varargin{:});
                
                % Now add your own code, such as domain constraints.
                % Exponential does not need any domain constraint.
                
                % Let YALMIP know about convexity etc
                varargout{1} = F;
                varargout{2} = struct('convexity','convex','monotonicity','increasing','definiteness','positive');
                varargout{3} = X;
                
            case 'milp'
                    varargout{1} = [];
                    varargout{2} = [];
                    varargout{3} = [];                
            otherwise
                error('SDPVAR/EXP called with CHAR argument?');
        end
    otherwise
        error('SDPVAR/EXP called with CHAR argument?');
end
