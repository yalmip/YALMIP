function varargout=fix(varargin)
%FIX (overloaded)

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
                
                % I cannot come up with a simple model for FIX that doesn't
                % introduce a binary to model negative and positivity of
                % the argument.
                d1 = binvar(1,1);
                d2 = binvar(1,1);
                [M,m] = derivebounds(X);
                F = [X>=m*(1-d1), X<=M*(1-d2),X-d1+0.00001 <=  t <= X+d2-0.0001,integer(t), d1+d2 == 1];
                
                varargout{1} = F;
                varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
                varargout{3} = X;
                
            otherwise
                error('SDPVAR/FIX called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/FIX');
end
