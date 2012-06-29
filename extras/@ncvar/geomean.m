function varargout = geomean(varargin)
%GEOMEAN (overloaded)
%
% t = GEOMEAN(X)
%
% For Hermitian matrix X, returns det(X)^(1/length(X))
%
% For real vector X, returns prod(X)^(1/length(X))
%
% This concave function is monotonically growing in det(P)
% for P>0, so  it can be used for maximizing det(P),
% or to add lower bound constraints on the determinant.
%
% When GEOMEAN is used in a problem, the domain constraint
% set(X>0) is automatically added to the problem.
%
% If you only use GEOMEAN as the objective function, you
% may want to consider GEOMEAN2 which may results in a
% slightly smaller SDP model.
%
% See also SDPVAR, GEOMEAN2, SUMK, SUMABSK

% Author Johan Löfberg
% $Id: geomean.m,v 1.1 2006-08-10 18:00:20 joloef Exp $

switch class(varargin{1})

    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        X = varargin{1};
        [n,m] = size(X);
        if is(varargin{1},'hermitian') | min(n,m)==1
            varargout{1} = yalmip('addextendedvariable',mfilename,varargin{:});
        else
            % Create one variable for each column            
            y = [];
            for i = 1:m
                index = (1+n*(i-1)):i*n;
                x = extsubsref(X,index);
                y = [y yalmip('addextendedvariable',mfilename,x)];
            end
            varargout{1} = y;
        end

    case 'char' % YALMIP send 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            t = varargin{2}; % Second arg is the extended operator variable
            X = varargin{3}; % Third arg and above are the args user used when defining t.

            % Extend if not power of 2
            n = length(X);
            m=2^ceil(log2(n));
            X0=X;
            F = set([]);
            if n ~= m
                d=m-n;
                if size(X,1)==size(X,2)
                    % model determinant
                    % Convert to a real problem
                    if is(X,'complex')
                        X = [real(X) imag(X);-imag(X) real(X)];
                        n = length(X);
                        % We will now get (detX)^2, so we need to pad
                        % differently to get (detX)^(1/n)
                        d=2*m-n;
                    end

                    D = tril(sdpvar(n,n));
                    delta = diag(D);
                    F = set([X D;D' diag(delta)] > 0);
                    p = 2^ceil(log2(n));
                    if 1
                        x = [delta;ones(d,1)*t];
                    else
                        % More efficient in detset, but not tested 
                        % sufficiently yet
                        x = delta;
                    end

                elseif size(X,1)>size(X,2)
                    x = [X;ones(d,1)*t];

                else 
                    x = [X ones(1,d)*t];
                end
            else
                if size(X,1)==size(X,2)
                    % model determinant
                    % Convert to a real problem
                    if is(X,'complex')
                        X = [real(X) imag(X);-imag(X) real(X)];
                        n = length(X);
                        % We will now get (detX)^2, so we need to pad
                        % differently to get (detX)^(1/n)
                        d=2*m-n;
                    end

                    D = tril(sdpvar(n,n));
                    delta = diag(D);
                    F = set([X D;D' diag(delta)] > 0);
                    p = 2^ceil(log2(n));
                    x = delta;
                   
                elseif size(X,1)>size(X,2)
                    x = X;

                elseif size(X,1)<size(X,2)
                    x = X;
                end
            end

            varargout{1} = F + detset(t,x);

            if (min(size(X))>1) & ishermitian(X0)
                varargout{2} = struct('convexity','concave','monotonicity','none','definiteness','positive');
            else
                varargout{2} = struct('convexity','concave','monotonicity','increasing','definiteness','positive');
            end
            varargout{3} = X0; %We have altered the original input, lucky we saved it!
        elseif isequal(varargin{1},'milp')

            varargout{1} = [];
            varargout{2} = [];
            varargout{3} = [];

        else
        end
    otherwise
end
