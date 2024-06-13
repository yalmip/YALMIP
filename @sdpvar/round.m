function varargout=round(varargin)
%ROUND (overloaded)

switch class(varargin{1})
    
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        
        x = varargin{1};
        dim = size(x);
        x = reshape(x,prod(dim),1);
        y = [];
        for i = 1:prod(dim)
            y = [y;yalmip('define',mfilename,extsubsref(x,i))];
        end
        y = reshape(y,dim);
        varargout{1} = y;
        
    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
        switch varargin{1}
            case {'milp','graph','exact'}
                % Description using epigraphs
                t = varargin{2};
                X = varargin{3};
                
                F = [X-0.5 <= t <= X+0.5, integer(t)];
                [M,m] = derivebounds(X);
                if ~isinf(m)
                    F = [F, round(m) <= t];
                end
                if ~isinf(M)
                    F = [F, t <= round(M)];
                end
                                 
                varargout{1} = F;
                varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
                varargout{3} = X;
                
            otherwise
                error('SDPVAR/ROUND called with CHAR argument?');
        end
    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end
