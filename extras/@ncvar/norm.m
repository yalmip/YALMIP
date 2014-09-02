function varargout = norm(varargin)
%NORM (overloaded)
%
% t = NORM(x,P)
%
% The variable t can only be used in convexity preserving
% operations such as t<0, max(t,y)<1, minimize t etc.
%
%    For matrices...
%      NORM(X) models the largest singular value of X, max(svd(X)).
%      NORM(X,2) is the same as NORM(X).
%      NORM(X,1) models the 1-norm of X, the largest column sum, max(sum(abs(X))).
%      NORM(X,inf) models the infinity norm of X, the largest row sum, max(sum(abs(X'))).
%      NORM(X,'fro') models the Frobenius norm, sqrt(sum(diag(X'*X))).
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

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        % SHOULD NEVER HAPPEN, THIS SHOULD BE CAUGHT BY BUILT-IN
        error('Overloaded SDPVAR/NORM CALLED WITH DOUBLE. Report error')

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        if nargin == 1
            varargout{1} = yalmip('addextendedvariable',mfilename,varargin{1},2);
        else
            switch varargin{2}
                case {1,2,inf,'inf','fro'}
                    varargout{1} = yalmip('addextendedvariable',mfilename,varargin{:});
                otherwise
                    error('norm(x,P) only supported for P = 1, 2, inf and ''fro''');
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

                % FIX : Exploit symmetry to create smaller problem
                switch p
                    case 1
                        z = sdpvar(size(X,1),size(X,2),'full');
                        if min(size(X))>1
                            if isreal(X)
                                F = (-z <= X <= z);
                            else
                                F = ([]);
                                for i = 1:size(X,1)
                                    for j = 1:size(X,2)
                                        xi = extsubsref(X,i,j);
                                        zi = extsubsref(z,i,j);
                                        F = F + (cone([real(xi);imag(xi)],zi));
                                    end
                                end
                            end
                            F = F + (sum(z,1) <= t);
                        else
                            if isreal(X)
                                F = (-z <= X <= z) + (sum(z) <= t);
                                [M,m] = derivebounds(X);
                                bounds(z,0,max(abs([M -m]),[],2));
                                bounds(t,0,sum(max(abs([M -m]),[],2)));
                            else
                                F = ([]);
                                for i = 1:length(X)
                                    xi = extsubsref(X,i);
                                    zi = extsubsref(z,i);
                                    F = F + (cone([real(xi);imag(xi)],zi));
                                end
                                F = F + (sum(z) <= t);
                            end
                        end
                    case 2
                        z = sdpvar(size(X,1),size(X,2));
                        if min(size(X))>1
                            F = ([t*eye(size(X,1)) X;X' t*eye(size(X,2))])>=0;
                        else
                            F = (cone(X(:),t));
                        end
                    case {inf,'inf'}
                        if min(size(X))>1
                            z = sdpvar(size(X,1),size(X,2),'full');
                            if isreal(X)
                                F = (-z <= X <= z);
                            else
                                F = ([]);
                                for i = 1:size(X,1)
                                    for j = 1:size(X,2)
                                        xi = extsubsref(X,i,j);
                                        zi = extsubsref(z,i,j);
                                        F = F + (cone([real(xi);imag(xi)],zi));
                                    end
                                end
                            end
                            F = F + (sum(z,2) <= t);
                        else
                            if isreal(X)
                                F = (-t <= X <= t);
                                [M,m,infbound] = derivebounds(X);
                                if ~infbound
                                    F = F + (0<=t<=max(M));
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
                    otherwise
                end
                varargout{1} = F;
                varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive');
                varargout{3} = X;
            case 'milp'

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
                         F = ([]);
                         positive = find(m >= 0);
                         negative = find(M <= 0);
                         
                         % d(find(positive)) = 1;
                         % d(find(negative)) = 0;                         
                         if ~isempty(positive)                             
                             d = subsasgn(d,struct('type','()','subs',{{positive}}),1);
                         end
                         if ~isempty(negative)                             
                             d = subsasgn(d,struct('type','()','subs',{{negative}}),0);
                         end
                         
                         F = F + (X <= M.*d)     + (2*m.*d     <= absX+X <= 2*M.*d);
                         F = F + (X >= m.*(1-d)) + (2*m.*(1-d) <= absX-X <= 2*M.*(1-d));
                         F = F + (t - sum(absX) == 0);
                         
                    else

                        if 0

                            %2^n cases
                            %e.g in 2d,
                            %norm([x;y],inf) = y, y>0, y>x.y>-x
                            %                = x, x>0, x>y.x>-y

                            n = length(X);
                            X     = reshape(X,n,1);
                            absX  = sdpvar(n,1);
                            d     = binvar(n,1);
                            [M,m] = derivebounds(X);

                            F = (sum(d)==0);

                            top = 1;
                            for i = 1:n
                                xi = extsubsref(X,i);
                                y = extsubsref(X,setdiff(1:n,i));
                                for sign_abs_largest_variable = -1:2:1
                                    di = extsubsref(d,top);
                                    for j = setdiff(1:n,i)
                                        y = extsubsref(X,j);
                                        for sign_other = -1:2:1
                                            F = F + (xi*sign_abs_largest_variable >= sign_other*y);
                                            F = F + (-M*100*(1-di) <= xi*sign_abs_largest_variable-t <= t+M*100*(1-di));
                                        end
                                    end
                                    top = top + 1;
                                end
                            end




                        else
                            % OLD
                            n = length(X);
                            X     = reshape(X,n,1);
                            absX  = sdpvar(n,1);
                            d     = binvar(n,1);
                            [M,m] = derivebounds(X);
                            F = ([]);
                            F = F + (X <= M.*d)     + (2*m.*d     <= absX+X <= 2*M.*d);
                            F = F + (X >= m.*(1-d)) + (2*m.*(1-d) <= absX-X <= 2*M.*(1-d));
                            M = max(M,-m);
                            d = binvar(n,1);
                            F = F + (sum(d)==1);
                            F = F + (absX <= t <= absX + M.*(1-d));

                            kk = [];
                            ii = [];
                            for i = 1:n
                                k = [1:1:i-1 i+1:1:n]';
                                ii = [ii;repmat(i,n-1,1)];
                                kk = [kk;k];
                                Mm = M(k);
                            end
                            xii = extsubsref(absX,ii);
                            dii = extsubsref(d,ii);
                            xkk = extsubsref(absX,kk);
                            F = F + (xkk <= xii+(M(kk)-m(ii)).*(1-dii));
                        end

                        %   for i = 1:n
                        %       xi = extsubsref(absX,i);
                        %       di = extsubsref(d,i);
                        %       for k = [1:1:i-1 i+1:1:n]
                        %           xk = extsubsref(absX,k);
                        %           F = F + (xk <= xi+M(k)*(1-di));
                        %       end
                        %   end

                    end
                    varargout{1} = F;
                    varargout{2} = struct('convexity','milp','monotonicity','milp','definiteness','positive');
                    varargout{3} = X;
                end
            otherwise
                error('SDPVAR/NORM called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/NORM');
end