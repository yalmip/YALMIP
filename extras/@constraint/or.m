function varargout = or(varargin)
%OR (overloaded)

% Prune empty clauses
varargin = {varargin{find(~cellfun('isempty',varargin))}};
% Models OR using a nonlinear operator definition
switch class(varargin{1})
    case 'char'
        z = varargin{2};
        X = varargin{3};
        Y = varargin{4};        
      
        F = ([]);                
        switch class(X)
            case 'sdpvar'
                x = X;
                xvars = getvariables(x);               
                allextvars = yalmip('extvariables');
                if (length(xvars)==1) & ismembcYALMIP(xvars,allextvars)
                    [x,F] = expandor(x,allextvars,F);
                end

            case 'constraint'
                x = binvar(1,1); 
                F = F + (implies_internal(x,X));
            otherwise
        end
        switch class(Y)
            case 'sdpvar'
                y = Y;
                yvars = getvariables(y);
                allextvars = yalmip('extvariables');
                if (length(yvars)==1) & ismembcYALMIP(yvars,allextvars)
                    [y,F] = expandor(y,allextvars,F);
                end
            case {'constraint','lmi'}
                y = binvar(1,1);
                F = F + (implies_internal(y,Y));
            otherwise
        end

        xy = [x y];
        [M,m] = derivebounds(z);
        if m>=1
            varargout{1} = F + (sum(xy) >= 1);
        else
            varargout{1} = F + (sum(xy) >= z) + (z >= xy) +(binary(z));
        end
        X = sdpvar(X);
        Y = sdpvar(Y);
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = [xy(:);X(:);Y(:)];

    case {'sdpvar','constraint'}             
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
            if ismembcYALMIP(xvars,allextvars)
                [x1,F] = expandor(x1,allextvars,F);
            end

        case 'constraint'
            x1 = binvar(1,1);
            F = F + (iff_internal(X1,x1));
        otherwise
    end
    switch class(X2)
        case 'sdpvar'
            x2 = X2;
            yvars = getvariables(x2);
            if ismembcYALMIP(yvars,allextvars)
                [x2,F] = expandor(x2,allextvars,F);
            end
        case 'constraint'
            x2 = binvar(1,1);
            F = F + (iff_internal(X2,x2));
        otherwise
    end
    x = [x1 x2];
end