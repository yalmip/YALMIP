function varargout = norm(varargin)
%NORM (overloaded)
%
% t = NORM(x,P)
%
% The variable t can only be used in convexity preserving
% operations such as t<=1, max(t,y)<=1, minimize t etc.
%
%    For matrices...
%      NORM(X)       models the largest singular value of X, max(svd(X)).
%      NORM(X,2)     is the same as NORM(X).
%      NORM(X,1)     models the 1-norm of X, the largest column sum, max(sum(abs(X))).
%      NORM(X,inf)   models the infinity norm of X, the largest row sum, max(sum(abs(X'))).
%      NORM(X,'inf') same as above
%      NORM(X,'fro') models the Frobenius norm, sqrt(sum(diag(X'*X))).
%      NORM(X,'nuc') models the Nuclear norm, sum of singular values.
%      NORM(X,'*')   same as above
%      NORM(X,'tv')  models the (isotropic) total variation semi-norm 
%    For vectors...
%      NORM(V) = norm(V,2) = standard Euclidean norm.
%      NORM(V,inf) = max(abs(V)).
%      NORM(V,1) = sum(abs(V))
%
% SEE ALSO SUMK, SUMABSK

%% ***************************************************
% This file defines a nonlinear operator for YALMIP
%
% It can take three different inputs
% For DOUBLE inputs, it returns standard double values
% For SDPVAR inputs, it generates an internal variable
%
% When first input is 'model' it returns the graph
% in the first output and structure describing some
% properties of the operator.

%% ***************************************************
switch class(varargin{1})
    
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        if nargin == 1
            varargout{1} = yalmip('define',mfilename,varargin{1},2);
        else
            switch varargin{2}
                case {1,2,inf,'inf','fro'}
                    varargout{1} = yalmip('define',mfilename,varargin{:});
                case 'tv'
                    if ~isreal(varargin{1})
                        error('Total variation norm not yet implemented for complex arguments');
                    end
                    if min(varargin{1}.dim)==1
                        varargout{1} = norm(diff(varargin{1}),1);
                        return
                    end
                    varargout{1} = yalmip('define','norm_tv',varargin{:});
                case {'nuclear','*'}
                    if min(size(varargin{1}))==1
                        varargout{1} = norm(varargin{1},1);
                    else
                        varargout{1} = yalmip('define','norm_nuclear',varargin{:});
                    end
                otherwise
                    if isreal(varargin{1}) & min(size(varargin{1}))==1 & isnumeric(varargin{2})
                        varargout{1} = pnorm(varargin{:});
                    else
                        error('norm(x,P) only supported for P = 1, 2, inf, ''fro'' and ''nuclear''');
                    end
            end
        end
        
    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        switch varargin{1}
            case 'graph'
                t = varargin{2};
                X = varargin{3};
                p = varargin{4};
                
                % Code below complicated by two things
                % 1: Absolute value for complex data -> cone constraints on
                %    elements
                % 2: SUBSREF does not call SDPVAR subsref -> use extsubsref.m
                                
                switch p
                    case 1
                        if issymmetric(X)
                            Z = sdpvar(size(X,1),size(X,2));
                        else
                            Z = sdpvar(size(X,1),size(X,2),'full');
                        end
                        if min(size(X))>1
                            if isreal(X)
                                z = reshape(Z,[],1);
                                x = reshape(X,[],1);
                                F = (-z <= x <= z);
                            else
                                F = ([]);
                                zvec = reshape(Z,1,[]);
                                xrevec=reshape(real(X),1,[]);
                                ximvec=reshape(imag(X),1,[]);
                                F = [F,cone([zvec;xrevec;ximvec])];
                            end
                            F = F + (sum(Z,1) <= t);
                        else
                            if isreal(X)
                                % Standard definition
                                % F = (-t <= X <= t);
                                X = reshape(X,[],1);
                                Z = reshape(Z,[],1);
                                Xbase = getbase(X);
                                Constant = find(~any(Xbase(:,2:end),2));
                                if ~isempty(Constant)
                                    % Exploit elements without any
                                    % decision variables
                                    r1 = ones(length(Z),1);
                                    r2 = zeros(length(Z),1);
                                    r1(Constant) = 0;
                                    r2(Constant) = abs(Xbase(Constant,1));
                                    Z = Z.*r1 + r2;
                                end
                                F = (-Z <= X <= Z) + (sum(Z) <= t);                                                               
                            else                                                                
                                F = (cone([reshape(Z,1,[]);real(reshape(X,1,[]));imag(reshape(X,1,[]))]));                                                              
                                F = F + (sum(Z) <= t);
                            end
                        end
                    case 2                    
                        if min(size(X))>1
                            F = ([t*eye(size(X,1)) X;X' t*eye(size(X,2))])>=0;
                        else
                            F = cone(X(:),t);
                        end
                    case {inf,'inf'}
                        if min(size(X))>1
                            Z = sdpvar(size(X,1),size(X,2),'full');
                            if isreal(X)
                                F = (-Z <= X <= Z);
                            else
                                F = ([]);
                                for i = 1:size(X,1)
                                    for j = 1:size(X,2)
                                        xi = extsubsref(X,i,j);
                                        zi = extsubsref(Z,i,j);
                                        F = F + (cone([real(xi);imag(xi)],zi));
                                    end
                                end
                            end
                            F = F + (sum(Z,2) <= t);
                        else
                            if isreal(X)
                                F = (-t <= X <= t);
                                [M,m,infbound] = derivebounds(X);
                                if ~infbound
                                    F = F + (0 <= t <= max(max(abs([m M]))));
                                end
                            else
                                F = ([]);
                                for i = 1:length(X)
                                    xi = extsubsref(X,i);
                                    F = F + (cone([real(xi);imag(xi)],t));
                                end
                            end
                        end
                    case 'fro'
                        X.dim(1)=X.dim(1)*X.dim(2);
                        X.dim(2)=1;
                        F = (cone(X,t));
                    case 'nuclear'
                        U = sdpvar(X.dim(2));
                        V = sdpvar(X.dim(1));
                        F = [trace(U)+trace(V) <= 2*t, [U X';X V]>=0];
                    case 'tv'                        
                        Dx = [diff(X,1,1);zeros(1,X.dim(2))];
                        Dy = [diff(X,1,2) zeros(X.dim(1),1)];
                        T = sdpvar(X.dim(1),X.dim(2),'full');
                        F = cone([reshape(T,1,[]);reshape(Dx,1,[]);reshape(Dy,1,[])]);
                        F = [F, sum(sum(T)) <= t];
                    otherwise
                end
                varargout{1} = F;
                varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive','model','graph');
                varargout{3} = X;
                
            case 'exact'
                
                t = varargin{2};
                X = varargin{3};
                p = varargin{4};
                if ~isreal(X) | isequal(p,2) | isequal(p,'fro') | min(size(X))>1 % Complex valued data, matrices and 2-norm not supported
                    varargout{1} = [];
                    varargout{2} = [];
                    varargout{3} = [];
                else
                    if p==1
                        X     = reshape(X,length(X),1);
                        absX  = sdpvar(length(X),1);
                        d     = binvar(length(X),1);
                        [M,m] = derivebounds(X);
                        
                        if all(abs(sign(m)-sign(M))<=1)
                            
                            % Silly convex case. Some coding to care of the
                            % case sign(0)=0...
                            d = ones(length(X),1);
                            d(m<0)=-1;
                            F = (t - sum(absX) == 0) + (absX == d.*X);
                        else
                            
                            F = ([]);
                            
                            
                            % Some fixes to remove trivial constraints
                            % which caused problems in a user m
                            positive = find(m >= 0);
                            negative = find(M <= 0);
                            fixed = find(m==M);
                            if ~isempty(fixed)
                                positive = setdiff(positive,fixed);
                                negative = setdiff(negative,fixed);
                            end
                            if ~isempty(positive)
                                d = subsasgn(d,struct('type','()','subs',{{positive}}),1);
                            end
                            if ~isempty(negative)
                                d = subsasgn(d,struct('type','()','subs',{{negative}}),0);
                            end
                            if ~isempty(fixed)
                                notfixed = setdiff(1:length(m),fixed);
                                addsum = sum(abs(m(fixed)));
                                m = m(notfixed);
                                M = M(notfixed);
                                X = extsubsref(X,notfixed);
                                absX = extsubsref(absX,notfixed);
                                d = extsubsref(d,notfixed);
                            else
                                addsum = 0;
                            end
                            
                            maxABSX = max([abs(m) abs(M)],[],2);
                            % d==0  ---> X<0 and absX = -X
                            F = F + (X <= M.*d)     + (0 <= absX+X <= 2*maxABSX.*d);
                            % d==1  ---> X>0 and absX = X
                            F = F + (X >= m.*(1-d)) + (0 <= absX-X <= 2*maxABSX.*(1-d));
                            
                            F = F + (t - sum(absX)-addsum == 0);
                        end
                        
                    else                        
                        F = max_integer_model([X;-X],t);                                                                       
                    end
                    varargout{1} = F;
                    varargout{2} = struct('convexity','convex','monotonicity','milp','definiteness','positive','model','integer');
                    varargout{3} = X;
                end
            otherwise
                error('SDPVAR/NORM called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/NORM');
end

function F = findmax(F,M,m,X,t)

n = length(X);
d = binvar(n,1);
F = F + (sum(d)==1);
F = F + (-(max(M)-min(m))*(1-d) <= t-X <= (max(M)-min(m))*(1-d));
kk = [];
ii = [];
for i = 1:n
    k = [1:1:i-1 i+1:1:n]';
    ii = [ii;repmat(i,n-1,1)];
    kk = [kk;k];
    Mm = M(k)-m(i);
end
xii = extsubsref(X,ii);
dii = extsubsref(d,ii);
xkk = extsubsref(X,kk);
F = F + (xkk <= xii+(M(kk)-m(ii)).*(1-dii));