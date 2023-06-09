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
        x = [];
        q = [];
        for j = 3:nargin
            X = varargin{j};
            switch class(X)
                case 'sdpvar'
                    xi = X;
                    q = [q;X(:)];
                    xvars = getvariables(x);               
                    allextvars = yalmip('extvariables');
                    if (length(xvars)==1) & ismembcYALMIP(xvars,allextvars)
                        [xi,F] = expandor(xi,allextvars,F);
                    end
                    x = [x xi];

                    case 'constraint'
                        xi = binvar(1,1); 
                        F = F + (implies_internal(xi,X));
                        x = [x xi];
                        q = [q;reshape(sdpvar(X),[],1)];
                otherwise
            end       
        end
        
        [M,m] = derivebounds(z);
        if m>=1
            varargout{1} = F + (sum(x) >= 1);
        else
            varargout{1} = F + (sum(x) >= z) + (z >= x) +(binary(z));
        end
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = [x(:);q];

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