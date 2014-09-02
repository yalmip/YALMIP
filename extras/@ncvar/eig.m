function varargout=eig(varargin)
%ABS (overloaded)
%
% y = eig(x)
%
% The variable y represents the sorted eigenvalues of x.
% 
% EIG is implemented using a nonlinear operator framework,
% but can currently not be used in any optmization model.



% ***************************************************
% This file defines a nonlinear operator for YALMIP
%
% It can take three different inputs
% For double inputs, it returns standard double values
% For sdpvar inputs, it genreates a an internal variable
% When first input is 'model' it generates the epigraph
%
% % ***************************************************
switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        % SHOULD NEVER HAPPEN, THIS SHOULD BE CAUGHT BY BUILT-IN
        error('Overloaded SDPVAR/EIG CALLED WITH DOUBLE. Report error')

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        X = varargin{1};
        if size(X,1)~=size(X,2)
            error('Matrix must be square.')
        end
        y=[];
        for i = 1:length(X)
            y=[y;yalmip('addextendedvariable','yeig',X,i)]; %ith eigenvalue
        end
        varargout{1}=y;            

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
        error('eig operator not supported');
        if isequal(varargin{1},'graph')
            t = varargin{2};
            X = varargin{3};
            varargout{1} = (-t <= X <= t);
            varargout{2} = 1; % Convex operator
            error('SDPVAR/ABS called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/ABS');
end
