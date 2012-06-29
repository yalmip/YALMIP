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
% $Id: sort.m,v 1.1 2006-08-10 18:00:22 joloef Exp $

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/SORT CALLED WITH DOUBLE. Report error')

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.

        if nargin > 1 | min(size(varargin{1}))>1
            error('SDPVAR/SORT only supports simple 1-D sorting'),
        end

        x = varargin{1};
        data.D = binvar(length(x),length(x),'full');
        data.V = sdpvar(length(x),length(x),'full');
        y = [];

        for i = 1:length(x)
            data.i = i;
            data.isthisloc = 0;
            y = [y;yalmip('addextendedvariable',mfilename,x,data)];%i,P,V)];
        end
        loc = [];
        for i = 1:length(x)
            data.i = i;
            data.isthisloc = 1;
            loc = [loc;yalmip('addextendedvariable',mfilename,x,data)];
        end
        varargout{1} = y;
        varargout{2} = loc;

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
        switch varargin{1}
            case {'milp','graph'}
                % Description using epigraphs
                t = varargin{2};
                X = varargin{3};
                data = varargin{4};
                %i = varargin{4};
                %D = varargin{5};
                %V = varargin{6};

                % Call external to allow subsrefs in classs
                [F,vars] = sort_internal(t,X,data);%i,D,V);

                varargout{1} = F;
                varargout{2} = struct('convexity','milp','monotonicity','milp','definiteness','milp');
                varargout{3} = X;

                % Currently a hack, gneral feature coming soon...
                varargout{2}.models = vars;

            otherwise
                error('SDPVAR/SORT called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/SORT');
end
