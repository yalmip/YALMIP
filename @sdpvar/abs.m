function varargout=abs(varargin)
%ABS (overloaded)

switch class(varargin{1})
    case 'double'
        error('Overloaded SDPVAR/ABS CALLED WITH DOUBLE. Report error')

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        if isreal(varargin{1})
            varargout{1} = yalmip('define',mfilename,varargin{1});
            B = getbase(varargin{1});
            if all(all(fix(B)==B))
                CurrInt = yalmip('quantvariables');
                if all(ismember(getvariables(varargin{1}),CurrInt))
                    yalmip('setintvariables',[CurrInt getvariables(varargout{1})]);
                end
            end
        else
            % For complex args, abs(X) is defined [norm(X(i,j),2)] in MATLAB
            y = [];
            x = varargin{1};
            for i = 1:size(x,1)
                temp = [];
                for j = 1:size(x,2)
                    if isempty(temp)
                        temp = norm(extsubsref(x,i,j));
                    else
                        temp = [temp norm(extsubsref(x,i,j))];
                    end
                end
                if isempty(y)
                    y = temp;
                else
                    y = [y;temp];
                end
            end
            varargout{1} = y;
        end

    case 'char' % YALMIP send 'graph' when it wants the epigraph or hypograph
        switch varargin{1}
            case 'graph'
                % Description using epigraphs
                t = varargin{2};
                X = varargin{3};             
                varargout{1} = [1 -1;-1 -1]*[X;t] <= [0;0];                
                varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','graph');
                varargout{3} = X;

            case {'exact','integer','callback'}
                % Exact description using binary variables
                t = varargin{2};
                X = varargin{3};
                d = varargin{4};
                [M,m]=derivebounds(X);
                if m>=0
                    F = (t == X);
                elseif M<=0
                    F = (t == -X);
                else
                    maxABSX = max([abs(m) abs(M)],[],2);
                    F = [[0 1 0;
                     0 -1 0;   
                     1 0 -M;
                     1 1 -2*maxABSX;
                     -1 -1 0;
                     -1 0 -m;
                     -1 1 2*maxABSX;1 -1 0]*[X;t;d] <= [maxABSX;0;0;0;0;-m;2*maxABSX;0]];
                end

                varargout{1} = F;
                varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','integer');
                varargout{3} = X;

            otherwise
                error('SDPVAR/ABS called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/ABS');
end