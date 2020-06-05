function varargout = SIGN(varargin)

switch class(varargin{1})

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
                if isequal(getbase(X),[0 1]) && ismember(getvariables(X),yalmip('tempintvariables'))
                    % Numerically unstable case with integer X
                    % d1: Negative <= -1
                    % d2: 0
                    % d3: Postive >= 1
                    F = [X >= m*d1, X <= M*d3, X <= (1-d1)*(M+1)-1,X >=(1-d3)*(m-1)+1, m*(1-d2) <= X <= M*(1-d2), t == -d1 + d3,d1+d2+d3==1];                    
                else
                    % Numerically unstable case with continuous X
                    F = [X >= m*d1, X <= M*d3, X <= (1-d1)*M,X >=(1-d3)*m, m*(1-d2) <= X <= M*(1-d2), t == -d1 + d3,d1+d2+d3==1];
                end

                varargout{1} = F;
                varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none');
                varargout{3} = X;

            otherwise
                error('SDPVAR/SIGN called with CHAR argument?');
        end
    otherwise
        error(['SDPVAR/' upper(mfilename) ' called with weird argument']);
end
