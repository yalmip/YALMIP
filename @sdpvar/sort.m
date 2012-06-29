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
% $Id: sort.m,v 1.9 2007-08-02 19:17:36 joloef Exp $

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
            y = [y;yalmip('define',mfilename,x,data)];%i,P,V)];
        end
        loc = [];
        for i = 1:length(x)
            data.i = i;
            data.isthisloc = 1;
            loc = [loc;yalmip('define',mfilename,x,data)];
        end
        varargout{1} = y;
        varargout{2} = loc;

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph

        t = varargin{2};
        X = varargin{3};
        data = varargin{4};

        % Call external to allow subsrefs in classs        
        [F,vars] = sort_internal(t,X,data);

        varargout{1} = F;
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = X;

        % Inofficial way to model several nonlinear variables in
        % one call
        varargout{2}.models = vars;
    otherwise
        error('Strange type on first argument in SDPVAR/SORT');
end
