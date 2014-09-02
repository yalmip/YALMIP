function varargout = pwq_yalmip(varargin)
%PWQ_YALMIP Defines a piecewise quadratic function using data from MPT
%
%Only intended for internal use in YALMIP
%
% Currently only a container for PWQ functions. Can not be
% used in actual optmization problem.

switch class(varargin{1})

    case {'struct','cell'} % Should only be called internally

        if isa(varargin{2},'double')
            % Called from YALMIP to get double
            pwastruct = varargin{1};
            x = varargin{2};
            index = varargin{5};

            val = inf;
            for i = 1:length(pwastruct)
                [ii,jj] = isinside(pwastruct{i}.Pn,x);
                if ii
                    for k = 1:length(jj)
                        Q = pwastruct{i}.Ai{jj(k)};
                        if index>1 | min(size(pwastruct{i}.Bi{jj(k)}))>1
                            % FIX: Why?? Where is this feature used
                            val = min(val,x'*Q*x + reshape(pwastruct{i}.Bi{jj(k)}(index,:),1,[])*x+pwastruct{i}.Ci{jj(k)}(index));
                        else
                            val = min(val,x'*Q*x + reshape(pwastruct{i}.Bi{jj(k)},1,[])*x+pwastruct{i}.Ci{jj(k)});
                        end
                    end
                end
            end
            if isinf(val)
                val = nan;
            end
            varargout{1} = min(val);
            return
        end

        if nargin<3
            pwaclass = 'general'
        end

        if isa(varargin{1},'struct')
            varargin{1} = {varargin{1}};
        end

        % Put in standard format
        if ~isfield(varargin{1}{1},'Bi')
            if ~isfield(varargin{1}{1},'Fi')
                error('Wrong format on input to PWQ (requires Bi or Fi)');
            else
                for i = 1:length(varargin{1})
                    varargin{1}{1}.Ai = cell(1, varargin{1}{i}.Fi);
                    varargin{1}{i}.Bi = varargin{1}{i}.Fi
                    varargin{1}{i}.Ci = varargin{1}{i}.Gi
                end
            end
        end

        if ~isfield(varargin{1}{1},'Pfinal')
            error('Wrong format on input to PWQ (requires field Pn)');
        end

        % This will be a container for binary variables in the furture
        varargin{end+1} = binvar(length(varargin{1}{1}.Pn),1);

        % Create one variable for each row
        % Inefficient but the only way currently in YALMIP
        varargout{1} = [];
        varargin{end+1} = 1;

        for i = 1:length(varargin{1}{1}.Ci{1})
            varargin{end} = i;
            varargout{1} = [varargout{1};yalmip('define',mfilename,varargin{:})];
        end

    case 'char' % YALMIP sends 'model' when it wants the epigraph or hypograph


        switch varargin{1}

            case 'graph'

                varargout{1} = [];
                varargout{2} = struct('convexity','none','monotoncity','none','definiteness','none');
                varargout{3} = [];

            case {'integer','exact'}

                % FIX : Should create case for overlapping convex PWAs,
                % used in a nonconvex fashion...

                % Can only generate the first class of PWA functions
                t     = varargin{2};     % The YALMIP variables modelling this pwa
                pwq_struct = varargin{3};% MPT structure
                x     = varargin{4};     % Argument
                flag  = varargin{5};     % Type of PWA function
                d     = varargin{6};     % Binary for nonconvex cases
                index = varargin{7};     % Which row in Bix+Ci

                switch flag

                    case {'general','convex'}

                        if length(d)==1
                            % Don't introduce any binary variables when
                            % there is only one quadratic cost
                            F = ([]);
                            cost = x'*pwq_struct{1}.Ai{1}*x+pwq_struct{1}.Bi{1}*x+pwq_struct{1}.Ci{1};
                        else
                            n = length(x);
                            m = length(d);
                            z = sdpvar(n,m,'full');
                            F = (sum(d) == 1);
                            cost = 0;
                            [Mm,mm] = derivebounds(x);
                            try
                                [aux,mx,Mx] = bounding_box(pwq_struct{1}.Pfinal);
                            catch
                                Mx = 10000*ones(length(x),1);
                                mx = -10000*ones(length(x),1);
                            end
                            Mx = min([Mx Mm],[],2);
                            mx = max([mx mm],[],2);
                            for i = 1:m
                                %    bounds(z(:,i),mx,Mx);
                                F = F + (-(Mx-mx)*(1-d(i)) <= z(:,i)-x <= (Mx-mx)*(1-d(i)));
                                F = F + (mx*d(i) <= z(:,i) <= Mx*d(i));
                                [H,K] = double(pwq_struct{1}.Pn(i));
                                [M,m]  = derivebounds(H*x-K);
                                F = F + (H*x-K <= M*(1-d(i)));
                                cost = cost + z(:,i)'*pwq_struct{1}.Ai{i}*z(:,i)+pwq_struct{1}.Bi{i}(:)'*z(:,i)+d(i)*pwq_struct{1}.Ci{i};
                            end
                        end

                        varargout{1} = F;
                        varargout{2} = struct('convexity','none','monotoncity','none','definiteness','none','model','integer');
                        varargout{2}.replacer = cost;
                        varargout{3} = x;

                    case {'convexoverlapping'}

                        n = length(x);
                        cost = 0;
                        r = binvar(length(pwq_struct),1);
                        d = {};

                        % Derive bounds from polytopes
                        [aux,mx,Mx] = bounding_box(pwq_struct{1}.Pfinal);
                        for i = 2:length(pwq_struct)
                            [aux,L,U] = bounding_box(pwq_struct{i}.Pfinal);
                            Mx = max([Mx U],[],2);
                            mx = min([mx L],[],2);
                        end
                        bounds(x,mx,Mx);

                        % Some overlapping function is active
                        % (it will automatically be the smallest one)
                        F = (sum(r) == 1);
                        for j = 1:length(pwq_struct)

                            if length(pwq_struct{j}.Pn) == 1
                                d{j} = r(j);
                            else
                                d{j} = binvar(length(pwq_struct{j}.Pn),1);
                                F = F + (sum(d{j}) == r(j));
                            end
                            m = length(d{j});
                            z = sdpvar(n,m,'full');

                            for i = 1:m
                                bounds(z(:,i),mx,Mx);
                                F = F + (-(Mx-mx)*(1-d{j}(i)) <= z(:,i)-x <= (Mx-mx)*(1-d{j}(i)));
                                F = F + (mx*d{j}(i) <= z(:,i) <= Mx*d{j}(i));
                                [H,K] = double(pwq_struct{j}.Pn(i));
                                [M,m]  = derivebounds(H*x-K);
                                F = F + (H*x-K <= M*(1-d{j}(i)));
                                cost = cost + z(:,i)'*pwq_struct{j}.Ai{i}*z(:,i)+pwq_struct{j}.Bi{i}*z(:,i)+d{j}(i)*pwq_struct{j}.Ci{i};
                            end
                        end

                        varargout{1} = F;
                        varargout{2} = struct('convexity','none','monotoncity','none','definiteness','none','model','integer');
                        varargout{2}.replacer = cost;
                        varargout{3} = x;

                    otherwise

                        varargout{1} = [];
                        varargout{2} = struct('convexity','convex','monotoncity','none','definiteness','none');
                        varargout{3} = [];
                        return
                end

            otherwise
                error('PWA_YALMIP called with CHAR argument?');
        end

    otherwise
        error('Strange type on first argument in PWQ_YALMIP');
end