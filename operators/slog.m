function varargout = slog(varargin)
%ENTROPY
%
% y = SLOG(x)
%
% Computes/declares concave shifted logarithm log(1+x)
%
% Implemented as evalutation based nonlinear operator. Hence, the concavity
% of this function is exploited to perform convexity analysis and rigorous
% modelling. Implemented in order to avoid singularities in logarithm
% evaluation.

switch class(varargin{1})
    case 'double'
        x = varargin{1};       
        varargout{1} = log(abs(1+x)+sqrt(eps));
        
    case 'sdpvar'         
        varargout{1} = check_for_special_cases(varargin{:});
        if isempty(varargout{1})
            varargout{1} = InstantiateElementWise(mfilename,varargin{:});
        end

    case 'char'

        X = varargin{3};
        F = (X >= -1+eps);

        operator = struct('convexity','concave','monotonicity','increasing','definiteness','none','model','callback');
        operator.convexhull = @convexhull;
        operator.range = [-inf inf];
        operator.domain = [-1 inf];
        operator.derivative = @(x) (1./(abs(1+x)+sqrt(eps)));
        operator.inverse = @(x)(exp(x)-1);
        
        varargout{1} = F;
        varargout{2} = operator;
        varargout{3} = X;

    otherwise
        error('SDPVAR/SLOG called with CHAR argument?');
end


function f = check_for_special_cases(g);
f = [];
% Check for slog(1+a(x)/y)
vars = getvariables(g);
[mt,variabletype] = yalmip('monomtable');
% All signomials
if all(variabletype(vars)==4)
    % All xi/yi
    local_mt = mt(vars,:);
    for i = 1:size(local_mt,1)
        if nnz(local_mt(i,:) < 0) ~= 1
            return
        end
        if nnz(local_mt(i,:)) >2
            return
        end
        if ~all(ismember(local_mt(i,:),[-1 0 1]))
            return
        end           
    end
    % OK, everything is xi/yi    
    for i = 1:length(g)
         gi = extsubsref(g,i);       
         [~,common] = find(mt(getvariables(gi),:) < 0);
         y = recover(common(1));
         x = gi*y;
         if length(g)==1
             % Testing some display issues
             f =slogfrac([x;y]);
         else
            f = [f;slogfrac([x;y])];
         end
    end    
else
    return
end


function [Ax, Ay, b, K] = convexhull(xL,xU)
K = [];
if 1+xL <= 0
    fL = inf;
else
    fL = log(1+xL);
end
fU = log(1+xU);
dfL = 1/(1+xL);
dfU = 1/(1+xU);
xM = (xU + xL)/2;
fM = log(1+xM);
dfM = 1/(1+xM);

[Ax,Ay,b] = convexhullConcave(xL,xM,xU,fL,fM,fU,dfL,dfM,dfU);
remove = isinf(b) | isinf(Ax) | isnan(b);
if any(remove)
    remove = find(remove);
    Ax(remove)=[];
    b(remove)=[];
    Ay(remove)=[];
end