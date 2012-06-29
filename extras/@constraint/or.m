function varargout = or(varargin)
%OR (overloaded)

% Author Johan Löfberg 
% $Id: or.m,v 1.11 2007-08-02 19:33:16 joloef Exp $   

% Models OR using a nonlinear operator definition
switch class(varargin{1})
    case 'char'
        z = varargin{2};
        X = varargin{3};
        Y = varargin{4};        
      
        F = set([]);                
        switch class(X)
            case 'sdpvar'
                x = X;
                xvars = getvariables(x);               
                allextvars = yalmip('extvariables');
                if (length(xvars)==1) & ismembc(xvars,allextvars)
                    [x,F] = expandor(x,allextvars,F);
                end

            case 'constraint'
                x = binvar(1,1); 
                F = F + set(implies_internal(x,X));
            otherwise
        end
        switch class(Y)
            case 'sdpvar'
                y = Y;
                yvars = getvariables(y);
                allextvars = yalmip('extvariables');
                if (length(yvars)==1) & ismembc(yvars,allextvars)
                    [y,F] = expandor(y,allextvars,F);
                end
            case 'constraint'
                y = binvar(1,1);
                F = F + set(implies_internal(y,Y));
            otherwise
        end

        xy = [x y];
        [M,m] = derivebounds(z);
        if m>=1
            varargout{1} = F + set(sum(xy) >= 1);
        else
            varargout{1} = F + set(sum(xy) > z) + set(z > xy) +set(binary(z));
        end
        varargout{2} = struct('convexity','none','monotonicity','exact','definiteness','none','model','integer');
        varargout{3} = xy;

    case {'sdpvar','constraint'}
        x = varargin{1};
        y = varargin{2};
        varargout{1} = yalmip('define','or',varargin{:});

    otherwise
end


function [x,F] = expandor(x,allextvars,F)

xmodel = yalmip('extstruct',getvariables(x));

if isequal(xmodel.fcn,'or')
    X1 = xmodel.arg{1};
    X2 = xmodel.arg{2};

    switch class(X1)
        case 'sdpvar'
            x1 = X1;
            xvars = getvariables(x1);
            if ismembc(xvars,allextvars)
                [x1,F] = expandor(x1,allextvars,F);
            end

        case 'constraint'
            x1 = binvar(1,1);
            F = F + set(iff_internal(X1,x1));
        otherwise
    end
    switch class(X2)
        case 'sdpvar'
            x2 = X2;
            yvars = getvariables(x2);
            if ismembc(yvars,allextvars)
                [x2,F] = expandor(x2,allextvars,F);
            end
        case 'constraint'
            x2 = binvar(1,1);
            F = F + set(iff_internal(X2,x2));
        otherwise
    end
    x = [x1 x2];
end