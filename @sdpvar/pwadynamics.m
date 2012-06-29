function varargout=sort(varargin)
%SORT (overloaded)
%
% [t,loc] = sort(x)
%
% The variable t will be the sorted version of x.
%
% SORT is implemented in the nonlinear operator framework using a big-M
% model.


% Author Johan Löfberg
% $Id: pwadynamics.m,v 1.3 2007-08-02 19:17:36 joloef Exp $

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/SORT CALLED WITH DOUBLE. Report error')

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.

        if 0
            error('')
        end

        y = [];
        
        update = varargin{1};
        for j = 1:length(update)
           y = [y;yalmip('define',mfilename,varargin{:},j)];
        end
        varargout{1} = y;        

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
        switch varargin{1}
            case {'milp','graph'}
                % Description using epigraphs
                t = varargin{2};
                X = varargin{3};
               

                % Call external to allow subsrefs in classs
                [F,vars] = pwadynamics_internal(t,varargin{3:end-1});

                varargout{1} = F;
                varargout{2} = struct('convexity','milp','monotonicity','milp','definiteness','milp');
                varargout{3} = X;

                % Inofficial way to model several nonlinear variables in
                % one call
                varargout{2}.models = vars;

            otherwise
                error('SDPVAR/SORT called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/SORT');
end
