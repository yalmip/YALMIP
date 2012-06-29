function varargout = lmixor(varargin)
%XOR (overloaded)

% Author Johan Löfberg
% $Id: lmixor.m,v 1.3 2007-07-29 17:32:28 joloef Exp $

% Models XOR using a nonlinear operator definition
switch class(varargin{1})
    case 'char'
        z = varargin{2};

        allextvars = yalmip('extvariables');
        X = {};
        for i = 3:nargin
            Xtemp = expandxor(varargin{i},allextvars);
            for j = 1:length(Xtemp)
                X{end + 1} = Xtemp{j};
            end
        end

        F = set([]);
        x = binvar(length(X),1);

        for i = 1:length(X)
            F = F + set(implies_internal(extsubsref(x,i),X{i}));
        end

        varargout{1} = F + set(sum(x) == 1);
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','extra','marker','model','integer');
        varargout{3} = recover(depends(F));

    case {'lmi'}
        x = varargin{1};
        y = varargin{2};
        varargout{1} = set(yalmip('define','lmior',varargin{:}) == 1);

    otherwise
end

function x = expandxor(x,allextvars)

if length(getvariables(x))>1
    x = {x};
    return;
end

xmodel = yalmip('extstruct',getvariables(x));

if ~isempty(xmodel) & isequal(xmodel.fcn,'lmior')
    x1 = xmodel.arg{1};
    x2 = xmodel.arg{2};
    if  ismembc(getvariables(x1),allextvars)
        x1 = expandxor(x1,allextvars);
    else
        x1 = {x1};
    end
    if  ismembc(getvariables(x2),allextvars)
        x2 = expandxor(x2,allextvars);
    else
        x2 = {x2};
    end
    x = {x1{:},x2{:}};
else
    x = {x};
end