function varargout = norm(varargin)
%SIGN (overloaded)

% Author Johan Löfberg
% $Id: sign.m,v 1.1 2006-08-10 18:00:22 joloef Exp $


%% ***************************************************
% This file defines a nonlinear operator for YALMIP
%
% It can take three different inputs
% For DOUBLE inputs, it returns standard double values
% For SDPVAR inputs, it generates an internal variable
%
% When first input is 'model' it returns the graph
% in the first output and structure describing some
% properties of the operator.

%% ***************************************************
switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        % SHOULD NEVER HAPPEN, THIS SHOULD BE CAUGHT BY BUILT-IN
        error('Overloaded SDPVAR/NORM CALLED WITH DOUBLE. Report error')

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        if length(varargin{1}) == 1
            varargout{1} = yalmip('addextendedvariable',mfilename,varargin{1});
        else
            y = [];
            n = size(varargin{1},1);
            m = size(varargin{1},2);
            varargin{1} = reshape(varargin{1},n*m,1);
            for i = 1:prod(size(varargin{1}))
                inparg = extsubsref(varargin{1},i);
                if isa(inparg,'double')
                    y = sign(inparg);
                else
                    y = [y yalmip('addextendedvariable',mfilename,inparg)];
                end
            end
            y = reshape(y,n,m);
            varargout{1} = y;
        end

    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        switch varargin{1}
            case 'graph'
                t = varargin{2};
                X = varargin{3};

                % No linear model available

                varargout{1} = [];
                varargout{2} = [];
                varargout{3} = [];
            case 'milp'

                t = varargin{2};
                X = varargin{3};
                p = varargin{4};

                d = binvar(1,1);
                [M,m] = derivebounds(X);
                F = set([]);
                F = set(X > d*m) + set(-2*d+1 <= t <= 1);
                F = set(X < (1-d)*M) + set(-1     <= t <= -1 + 2*(1-d));

                varargout{1} = F;
                varargout{2} = struct('convexity','milp','monotonicity','milp','definiteness','positive');
                varargout{3} = X;

            otherwise
                error('SDPVAR/NORM called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/NORM');
end
