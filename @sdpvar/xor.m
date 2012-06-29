function varargout = xor(varargin)
%XOR (overloaded)
%
%    z = xor(x1,x2,...,xn)
%
% The OR operator is implemented using the concept of nonlinear operators
% in YALMIP. XOR defines a new so called derived variable that can be
% treated as any other variable in YALMIP. When SOLVESDP is issued,
% constraints are added to the problem to model the XOR operator. The new
% constraints add constraints to ensure that z and x satisfy the
% truth-table for XOR.
%
% It is assumed that x are binary variables (either explicitely
% declared using BINVAR, or constrained using BINARY.)
%
% See also SDPVAR/AND, SDPVAR/OR, BINVAR, BINARY

% Author Johan Löfberg
% $Id: xor.m,v 1.3 2007-07-29 17:32:29 joloef Exp $

% Models OR using a nonlinear operator definition
switch class(varargin{1})
    case 'char'
        z = varargin{2};

        allextvars = yalmip('extvariables');
        xy = [];
        for i = 3:nargin
            x{i-2} = varargin{i};
            xvars = getvariables(x{i-2});
            if (length(xvars)==1) & ismembc(xvars,allextvars)
                x = expandor(x,allextvars);
            end
            xy = [xy x{i-2}];
        end

        n = length(xy);

        [M,m] = derivebounds(z);

        % If user has constrained the XOR operator to true, we can add that
        % constraint very easily
        if m>=1
            varargout{1} = set(sum(xy) == 1);
        else
            T1 = -ones(n);
            T1 = T1 + 2*diag(ones(n,1));
            t = combnk(1:n,2);
            T2 = zeros(size(t,1),n);
            for i = 1:size(t,1)
                T2(i,t(i,1)) = 1;
                T2(i,t(i,2)) = -1;
            end
            varargout{1} = set(2 - T2*reshape(xy,n,1) >= z) + set(z >= T1*reshape(xy,n,1)) +set(binary(z));
        end

        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = xy;

    case 'sdpvar'
        x = varargin{1};
        y = varargin{2};
        varargout{1} = yalmip('addextendedvariable','xor',varargin{:});

    otherwise
end

function x = expandor(x,allextvars)

xmodel = yalmip('extstruct',getvariables(x));

if isequal(xmodel.fcn,'or')
    x1 = xmodel.arg{1};
    x2 = xmodel.arg{2};
    if  ismembc(getvariables(xmodel.arg{1}),allextvars)
        x1 = expandor(xmodel.arg{1},allextvars);
    end
    if  ismembc(getvariables(xmodel.arg{2}),allextvars)
        x2 = expandor(xmodel.arg{2},allextvars);
    end
    x = [x1 x2];
end