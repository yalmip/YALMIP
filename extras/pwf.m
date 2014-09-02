function varargout = pwf(varargin)
%PWF Defines a piecewise function
%
% t = PWF(h1,F1,h2,F2,...,hn,Fn)
%
% The variable t can only be used in convexity/concavity preserving
% operations, depending on the convexity of hi

switch class(varargin{1})
        
    case {'cell','struct'}
        % This should be an internal call to setup a PWA/PWQ function
        % defined in MPT 
                
        if nargin<3
            pwaclass = 'general'
        end

        if isa(varargin{1},'struct')
            varargin{1} = {varargin{1}};
        end

        % Put in standard format
        if ~isfield(varargin{1}{1},'Bi')
            if ~isfield(varargin{1}{1},'Fi')
                error('Wrong format on input to PWA (requires Bi or Fi)');
            else
                for i = 1:length(varargin{1})
                    varargin{1}{1}.Ai = cell(1, varargin{1}{i}.Fi);
                    varargin{1}{i}.Bi = varargin{1}{i}.Fi
                    varargin{1}{i}.Ci = varargin{1}{i}.Gi
                end
            end
        end
        
        if isempty(varargin{1}{1}.Ai{1})
            varargout{1} = pwa_yalmip(varargin{:});
        else
            if  nnz([varargin{1}{1}.Ai{:}]) == 0
                varargout{1} = pwa_yalmip(varargin{:});
            else
                varargout{1} = pwq_yalmip(varargin{:});
            end
        end
        
    case {'sdpvar','double'} % Overloaded operator for SDPVAR objects. Pass on args and save them.
        if length(varargin{1}) == 1
            
            % If it is a quadratic function, we might treat it more
            % efficiently by writing it as sum(d_i q_i(x)) with d_i binary.
            % This is implemented in the pwq_operator
            quadratic_objective = 1;
            linear_constraint =   1;
            for i = 1:nargin/2
                ok = 0;
                if isa(varargin{2*i-1},'sdpvar')
                    if isa(varargin{2*i},'constraint')
                        varargin{2*i} = lmi(varargin{2*i});
                    end
                    if isa(varargin{2},'lmi')
                        if is(varargin{2*i-1},'quadratic') |  is(varargin{2*i-1},'linear')
                            if all(is(varargin{2*i},'elementwise') | is(varargin{2*i},'elementwise'))
                                if is(sdpvar(varargin{2*i}),'linear')
                                    ok = 1;
                                end
                            end
                        end
                    end
                end
                if ~ok
                    break
                end
            end
            if ok
                z = [];
                for i = 1:nargin/2
                    z = [z depends(varargin{2*i-1}) depends(varargin{2*i})];
                end
                z = recover(unique(z));
                
                pwq{1}.Ai = {};
                pwq{1}.Bi = {};
                pwq{1}.Ci = {};
                pwq{1}.Pn = [];
                for i = 1:nargin/2
                    [Q,c,f,x,info] = quaddecomp(varargin{2*i-1},z);
                    pwq{1}.Ai = {pwq{1}.Ai{:},Q};
                    pwq{1}.Bi = {pwq{1}.Bi{:},c(:)};
                    pwq{1}.Ci = {pwq{1}.Ci{:},f};
                    F = getbase([sdpvar(varargin{2*i});sum(z)]);
                    pwq{1}.Pn = [pwq{1}.Pn polytope(-F(1:end-1,2:end),F(1:end-1,1))];                    
                end
                pwq{1}.Pfinal = union(pwq{1}.Pn);                                              
                varargout{1} = pwq_yalmip(pwq,z,'general');
            else
                varargout{1} = yalmip('define',mfilename,varargin{:});
            end
        else
            y = [];
            VV = varargin;
            for i = 1:length(varargin{1})
                VV = varargin;
                for j = 1:(length(varargin)/2)
                    VV{2*j-1} = VV{2*j-1}(i);
                end
                y = [y;yalmip('define',mfilename,VV{:})];
            end
            varargout{1} = y;
        end

    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph
        if isequal(varargin{1},'graph')
            varargout{1} = [];
            varargout{2} =  struct('convexity','failure','monotoncity','failure','definiteness','failure');;
            varargout{3} = [];
        elseif isequal(varargin{1},'milp')
            % pwf represented by t
            t = varargin{2};
            % Get functions and guards
            for i = 3:2:nargin
                f{(i-1)/2} = varargin{i};
                Guard{(i+1)/2-1} = varargin{i+1};
            end

            % Indicator for where we are
            indicators = binvar(length(f),1);
            % We are in some of the regions
            F = (sum(indicators) == 1);

            % Control where we are using the indicators and big-M
            X = [];
            for i = 1:length(Guard)
                Xi = sdpvar(Guard{i});
                if is(Xi,'linear')
                    [M,m] = derivebounds(Xi);
                else
                    m = -1e4;
                end
                F = F + (Xi >= m*(1-indicators(i)));
                X = [X;recover(depends(Guard{i}))];
            end
            % cost = sum(delta_i*cost_i(x)) = sum(cost_i(delta_i*x))
            X = recover(getvariables(X));

            if 0
                cost = 0;
                for i = 1:length(f)
                    [Q{i},c{i},g{i},xi,info] = quaddecomp(f{i},X);
                    if info
                        error('Only convex quadratic functions allowed in PWF');
                    end
                    [M,m] = derivebounds(xi);
                    z{i} = sdpvar(length(xi),1);
                    cost = cost + z{i}'*Q{i}*z{i}+c{i}'*z{i} + g{i}*indicators(i);
                    F = F + (xi+(1-indicators(i))*m<=z{i}<xi+(1-indicators(i)).*M);
                    F = F + (m*indicators(i) <= z{i} <= M*indicators(i));
                end
                F = F + (cost <= t);
            else
                fmax_max = -inf;
                fmin_min = inf;
                for i = 1:length(f)
                    [Q{i},c{i},g{i},xi,info] = quaddecomp(f{i},X);
                    if info
                        error('Only convex quadratic functions allowed in PWF');
                    end                    
                    [fmax,fmin] = derivebounds(c{i}'*xi+g{i});                                        
                    fmax_max = max(fmax,fmax_max);
                    fmin_min = min(fmin,fmin_min);
                    F = F + (c{i}'*xi +g{i} <= t + fmax*(1-indicators(i)));
                end
            end
            % To help variable strenghtening, constraint removal etc in
            % other functions, let us update the internal bound info
            bounds(t,fmin_min,fmax_max);

            % This function is convex when we relax the binary variables
            varargout{1} = F;
            varargout{2} = struct('convexity','milp','monotoncity','milp','definiteness','none');
            varargout{3} = X;
        else
            error('SDPVAR/NORM called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/NORM');
end