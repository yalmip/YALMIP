function varargout = sqrt(varargin)
%SQRT (overloaded)
%
% t = sqrt(x)
%
% The variable t can only be used in concavity preserving
% operations such as t>1, max t etc.
%
% When SQRT is used in a problem, the domain constraint
% set(x>0) is automatically added to the problem.
%
% See also SDPVAR, SDPVAR/GEOMEAN

% Author Johan Löfberg
% $Id: sqrt.m,v 1.15 2009-10-02 09:09:13 joloef Exp $

switch class(varargin{1})

    case 'double' % What is the numerical value of this argument (needed for displays etc)
        X = sqrt(varargin{1});

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.

        X = varargin{1};
        [n,m] = size(X);
        if is(varargin{1},'real') %& (n*m==1)
            varargout{1} = InstantiateElementWise(mfilename,varargin{:});
        else
            error('SQRT can only be applied to real scalars');
        end
      
        
    case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
        switch varargin{1}
            case 'graph'
            
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.
            if is(X,'linear')
                varargout{1} = set(cone([(X-1)/2;t],(X+1)/2));
                varargout{2} = struct('convexity','concave','monotonicity','increasing','definiteness','positive');
                varargout{3} = X;
            elseif is(X,'quadratic')
                [F,x] = check_for_special_cases(X,t);
                if isempty(F)
                   varargout{1} = [];
                varargout{2} = [];
                varargout{3} = [];
                else
                    varargout{1} = F;
                    varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','positive');    
                    varargout{3} = x;
                end
                
            else
                varargout{1} = [];
                varargout{2} = [];
                varargout{3} = [];
            end
            
            case 'callback'
                
                X = varargin{3};
                [U,L] = derivebounds(X);
                if L<0
                    F = set(X >= 0);
                else
                    F = [];
                end
                
                operator = struct('convexity','concave','monotonicity','increasing','definiteness','positive','model','callback');
                operator.range = [0 inf];
                operator.domain = [0 inf];
                operator.bounds = @bounds;
                operator.convexhull = @convexhull;
                operator.derivative = @(x) (1./(2*sqrt(abs(x)+eps)));%Yea, garbage for x<0
                
                varargout{1} = F;
                varargout{2} = operator;
                varargout{3} = X;
            otherwise
                varargout{1} = [];
                varargout{2} = [];
                varargout{3} = []; 
        end
    otherwise
end


function [F,x] = check_for_special_cases(q,t)
% Check if user is constructing sqrt(quadratic). If that is the case,
% return norm(linear)
F = [];
x = [];
if length(q)>1
    return
end
[Q,c,f,x,info] = quaddecomp(q);
if info==0 & nnz(Q)>0
    index = find(any(Q,2));
    if length(index) < length(Q)
        Qsub = Q(index,index);
        [Rsub,p]=chol(Qsub);
        if p==0
            [i,j,k] = find(Rsub);
            R = sparse((i),index(j),k,length(Qsub),length(Q));
        else
            R = [];
        end
    else
        [R,p]=chol(Q);
        if p & min(eig(full(Q)))>=-1e-12
            [U,S,V] = svd(full(Q));
            r = max(find(diag(S)));
            R = sqrtm(S(1:r,1:r))*V(:,1:r)';
            p = 0;
        end
    end
    d = 0.5*(R'\c);
    if p==0 &  f-d'*d>-1e-12
        F = cone([R*x+d;sqrt(f-d'*d)],t);
    end
end



function [L,U] = bounds(xL,xU)
if xL < 0
    % The variable is not bounded enough yet
    L = 0;
else
    L = sqrt(xL);
end
if xU < 0
    % This is an infeasible problem
    L = inf;
    U = -inf;
else
    U = sqrt(xU);
end

function [Ax, Ay, b] = convexhull(xL,xU)
if xL < 0 | xU == 0
    Ax = []
    Ay = [];
    b = [];
else
    xM = (1/(2*(sqrt(xU)-sqrt(xL))/(xU-xL)))^2;
    fL = sqrt(xL);
    fM = sqrt(xM);
    fU = sqrt(xU);
    dfL = 1/(2*sqrt(xL));
    dfM = 1/(2*sqrt(xM));
    dfU = 1/(2*sqrt(xU));
    [Ax,Ay,b] = convexhullConcave(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
    remove = isinf(b) | isinf(Ax) | isnan(b);
    if any(remove)
        remove = find(remove);
        Ax(remove)=[];
        b(remove)=[];
        Ay(remove)=[];
    end
end