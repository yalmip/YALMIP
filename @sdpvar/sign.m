function varargout = SIGN(varargin)
%SIGN (overloaded)

switch class(varargin{1})

    case 'double'
        error('Overloaded SDPVAR/SIGN CALLED WITH DOUBLE. Report error')

    case 'sdpvar' 
        varargout{1} = InstantiateElementWise(mfilename,varargin{:});

    case 'char' 
        switch varargin{1}
            case {'graph','exact'}

                t = varargin{2};
                X = varargin{3};
                
                d1 = binvar(1);
                d2 = binvar(1);
                d3 = binvar(1);
                [M,m] = derivebounds(X);
                %F = [X >= d*m,X <=(1-d)*M, t == 1-2*d];
                F = [X <= (1-d1)*M,X >=(1-d3)*m, m*(1-d2) <= X <= M*(1-d2), t == -d1 + d3,d1+d2+d3==1];

                varargout{1} = F;
                varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none');
                varargout{3} = X;

            otherwise
                error('SDPVAR/SIGN called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/SIGN');
end
