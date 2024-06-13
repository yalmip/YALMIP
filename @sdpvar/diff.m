function Y=diff(varargin)
%DIFF (overloaded)

X = varargin{1};
n = X.dim(1);
m = X.dim(2);

switch nargin
    case 1
        if n == 1
            % Default is diff along first dimension, unless we only have
            % one row. This happens to be the case in which the core function operates            
            Y = differ(X);
        else
            % MATLAB stadnard, just diff along first dimension. Reuse core
            % code by transposing matrix 
            Y = differ(X')';
        end            
        return

    case 2       
        % Note that MATLAB diff(randn(2,4),2) will cause a diff in two
        % different directions. 
        if isa(varargin{2},'sdpvar')
            error('I think you meant to use the JACOBIAN operator');
        end
        if n == 1
            if varargin{2}==1 || isempty(varargin{2})
                Y = differ(X);
            else
                % Recursive use
                Y = diff(differ(X),varargin{2}-1);
            end
        else
            if varargin{2}==1 || isempty(varargin{2})
                Y = differ(X')';
            else
                % Recursive use
                Y = diff(differ(X')',varargin{2}-1);
            end           
        end
        
    case 3
        if isa(varargin{2},'sdpvar')
            error('I think you meant to use the JACOBIAN operator');
        end
        if X.dim(varargin{3})-varargin{2} <=0
            dim = X.dim;
            dim(varargin{3}) = 0;
            Y = zeros(dim);
        elseif varargin{2}==1            
            if varargin{3} == 2
                Y = differ(X);
            elseif varargin{3} == 1 || isempty(varargin{3})
                Y = differ(X')';
            end
        elseif varargin{2} > 1
            % Recursive call for higher-order diff
            Y = diff(diff(X,varargin{2}-1,varargin{3}),1,varargin{3});
            return
        end
    otherwise
        error('To many input arguments.')
end


function Y = differ(X)

n = X.dim(1);
m = X.dim(2);
Y = X;
shift = [-speye(m-1) spalloc(m-1,1,0)] + [spalloc(m-1,1,0) speye(m-1)];
shift = kron(shift,speye(n));
Y.basis = shift*X.basis;
Y.dim(1) = n;
Y.dim(2) = m - 1;
% Reset info about conic terms
Y.conicinfo = [0 0];
Y = clean(Y);