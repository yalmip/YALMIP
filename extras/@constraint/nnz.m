function varargout = nnz(varargin)
%NNZ (overloaded)
%
%    n = nnz(X)
%
% NNZ applied to a constraint counts the number of satisfied constraints.
%
% The NNZ operator is implemented using the concept of nonlinear operators
% in YALMIP. NNZ(X) creates a new so called derived variable that can be
% treated as any other variable in YALMIP. When SOLVESDP is issued,
% logic constraints are added to the problem to model the NNZ operator.

switch class(varargin{1})

    case 'constraint'
        varargout{1} = yalmip('define','nnz',varargin{1});

    case 'char'
        t = varargin{2};
        x = varargin{3};
        switch varargin{1}
            case 'graph'
                x = lmi(x);
                if is(x,'elementwise')
                    x = sdpvar(x);
                    x = x(:);
                    [M,m] = derivebounds(x);
                    delta = binvar(length(x),1);
                    F = (x >= m.*(1-delta)) + (sum(delta) == t);
                elseif is(x,'equality')
                    x = sdpvar(x);
                    x = x(:);
                    [M,m,infbound] = derivebounds(x);
                    if infbound
                        warning('You have unbounded variables in IMPLIES leading to a lousy big-M relaxation.');
                    end
                    delta = binvar(length(x),1);
                    F = (M.*(1-delta) >= x >= m.*(1-delta)) + (sum(delta) == t);
                else
                    error('Constraint type not supported in cardinality. Make a feature request');
                end
                varargout{1} = F;
                varargout{2} = struct('convexity','concave','monotonicity','none','definiteness','none');
                varargout{3} = x;

            otherwise
                x = lmi(x);
                if is(x,'elementwise')
                    x = sdpvar(x);
                    x = x(:);
                    [M,m,infbound] = derivebounds(x);
                    if infbound
                            warning('You have unbounded variables in NNZ leading to a lousy big-M relaxation.');
                    end
                    n = length(x);
                    deltaD = binvar(n,1);
                    deltaU = binvar(n,1);
                    if  all(ismember(getvariables(x),[yalmip('intvariables') yalmip('binvariables')]))
                        % If an integer constraint not holds, it has to  be
                        % epsilon violated
                        F = ((M+1e-4).*(1-deltaU)-1e-4 >= x >= m.*(1-deltaD));
                    else
                        % User has to think harder him self
                        F = ((M).*(1-deltaU) >= x >= m.*(1-deltaD));
                    end
                    F = F + (sum(deltaD) == t);
                    F = F + (deltaU + deltaD == 1);
                    F = F + (0 <= t <= n);
                elseif is(x,'equality')
                    x = sdpvar(x);
                    x = x(:);
                    n = length(x);
                    [M,m] = derivebounds(x);
                    deltaZ   = binvar(n,1);
                    deltaNZU = binvar(n,1);
                    deltaNZD = binvar(n,1);
                    eps = 1e-4;
                    % deltaZ forces x to be zero
                    % deltaZU 1 is x > eps
                    % deltaZU 1 is x < -eps
                    F = (M.*(1-deltaZ) >= x >= m.*(1-deltaZ));
                    F = ((M+eps).*deltaNZU-eps>= x >= eps + (m-eps).*deltaNZD);
                    F = F + (sum(deltaZ) == t);
                    F = F + (deltaZ + deltaNZU + deltaNZD == 1);
                    F = F + (0 <= t <= n);
                else
                    error('Constraint type not supported in cardinality. Make a feature request');
                end
                varargout{1} = F;
                varargout{2} = struct('convexity','milp','monotonicity','milp','definiteness','milp');
                varargout{3} = x;
        end
    otherwise
end
