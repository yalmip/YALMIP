function varargout = norm_tv(varargin)
% NORM_TV Returns total variation semi-norm

switch class(varargin{1})
    
    case 'double' % What is the numerical value of this argument (needed for displays etc)
       
        X = varargin{1};
        [n,m] = size(X);
        Dx = [diff(X,1,1);zeros(1,m)];
        Dy = [diff(X,1,2) zeros(n,1)];
        Z = [Dx(:);Dy(:)];
        varargout{1} = sum(sqrt(sum(Z.^2,1)));                                                 
        
    case 'sdpvar'
        varargout{1} = yalmip('define',mfilename,varargin{:});
        
    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2};
            X = varargin{3};
            
            [n,m] = size(X);
            Dx = [diff(X,1,1);zeros(1,m)];
            Dy = [diff(X,1,2) zeros(n,1)];
            T = sdpvar(n,m,'full');
            F = cone([reshape(T,1,[]);reshape(Dx,1,[]);reshape(Dy,1,[])]);
            F = [F, sum(sum(T)) <= t];
            
            varargout{1} = F;
            varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','graph');
            varargout{3} = X;
        else
            varargout{1} = [];
            varargout{2} = [];
            varargout{3} = [];
        end
    otherwise
end