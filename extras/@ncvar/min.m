function y=min(varargin)
%MIN (overloaded)
%
% t = MIN(X,Y,DIM)
%
% Creates an internal structure relating the variable t with concave
% operator MIN(X).
%
% The variable t is primarily meant to be used in convexity preserving
% operations such as t>0, maximize t etc.
%
% If the variable is used in a non-convexity preserving operation, such as
% t<0, a mixed integer model will be derived.
%
% See built-in MIN for syntax.

% Author Johan Löfberg
% $Id: min.m,v 1.1 2006-08-10 18:00:21 joloef Exp $

% MIN is implemented as a nonlinear operator.
% However, for performance issues, it is not
% implemented in the default high-level way.
%
% The return of the double value and the
% construction of the epigraph/milp is done
% in the file model.m
%
% To study a better example of how to create
% your own nonlinear operator, check the
% function sdpvar/norm instead

% To simplify code flow, code for different #inputs
switch nargin
    case 1
        % Three cases:
        % 1. One scalar input, return same as output
        % 2. A vector input should give scalar output
        % 3. Matrix input returns vector output
        X = varargin{1};

        if max(size(X))==1
            y = X;
            return
        elseif min(size(X))==1
            y = yalmip('addextendedvariable','min',X);
            return
        else
            % This is just short hand for general command
            y = min(X,[],1);
        end

    case 2

        X = varargin{1};
        Y = varargin{2};
        [nx,mx] = size(X);
        [ny,my] = size(Y);
        if ~((nx*mx==1) | (ny*my==1))
            % No scalar, so they have to match
            if ~((nx==ny) & (mx==my))
                error('Array dimensions must match.');
            end
        end

        % Convert to compatible matrices
        if nx*mx==1
            X = X*ones(ny,my);
            nx = ny;
            mx = my;
        elseif ny==my
            Y = Y*ones(nx,mx);
            ny = nx;
            my = mx;
        end

        % Ok, done with error checks etc.
        % Introduce one extended variable/element
        % Lame code since we cannot call overloaded
        % subsref from inside of SDPVAR object...
        % This is most likely not a often used syntax,
        % so performance is neglected here...
        y = [];
        for i = 1:nx
            temp = [];
            for j = 1:mx
                Xij = extsubsref(X,i,j);
                Yij = extsubsref(Y,i,j);
                yij = min([Xij Yij]);
                temp = [temp yij];
            end
            y = [y;temp];
        end

    case 3

        X = varargin{1};
        Y = varargin{2};
        DIM = varargin{3};

        if ~(isa(X,'sdpvar') & isempty(Y))
            error('MIN with two matrices to compare and a working dimension is not supported.');
        end

        if ~isa(DIM,'double')
            error('Dimension argument must be 1 or 2.');
        end

        if ~(length(DIM)==1)
            error('Dimension argument must be 1 or 2.');
        end

        if ~(DIM==1 | DIM==2)
            error('Dimension argument must be 1 or 2.');
        end

        if DIM==1
            % Create one extended variable per column
            y = [];
            for i = 1:size(X,2)
                inparg = extsubsref(X,1:size(X,1),i);
                if isa(inparg,'sdpvar')
                    y = [y yalmip('addextendedvariable','min',inparg)];
                else
                    y = [y min(inparg)];
                end
            end
        else
            % Re-use code recursively
            y = min(X',[],1)';
        end

    otherwise
        error('Too many input arguments.')
end