function varargout=abs(varargin)
%ABS (overloaded)
%
% t = abs(x)
%
% The variable t can only be used in convexity preserving
% operations such as t<0, min t etc.

% Author Johan Löfberg 
% $Id: abs.m,v 1.1 2006-08-10 18:00:19 joloef Exp $

%% ***************************************************
% This file defines a nonlinear operator for YALMIP
%
% It can take three different inputs
% For double inputs, it returns standard double values
% For sdpvar inputs, it genreates a an internal variable
% When first input is 'model' it generates the epigraph

%% ***************************************************
switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/ABS CALLED WITH DOUBLE. Report error')

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        if isreal(varargin{1})
            varargout{1} = yalmip('addextendedvariable',mfilename,varargin{1});
        else
            % For complex args, abs(X) is defined [norm(X(i,j),2)] in MATLAB
            y = [];
            x = varargin{1};
            for i = 1:size(x,1)
                temp = [];
                for j = 1:size(x,2)
                    temp = [temp norm(extsubsref(x,i,j))];
                end
                y = [y;temp];
            end
            varargout{1} = y;
        end

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
        switch varargin{1}
            case 'graph'
                % Description using epigraphs
                t = varargin{2};
                X = varargin{3};
                varargout{1} = set(-t <= X <= t);
                varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive');
                varargout{3} = X;

            case 'milp'
                % Exact description using binary variables
                t = varargin{2};
                X = varargin{3};
                F = set([]);
                [M,m]=derivebounds(X);
                if m>=0
                    F = F + set(t == X);
                elseif M<0
                    F = F + set(t == -X);
                else
                    d = binvar(1,1);
                    F = F + set(X <= M*d)     + set(-2*(M-m)*d     <= t+X <= 2*(M-m)*d);
                    F = F + set(X >= m*(1-d)) + set(-2*(M-m)*(1-d) <= t-X <= 2*(M-m)*(1-d));
                end
                varargout{1} = F;
                varargout{2} = struct('convexity','milp','monotonicity','milp','definiteness','positive');
                varargout{3} = X;
            otherwise
                error('SDPVAR/ABS called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/ABS');
end
