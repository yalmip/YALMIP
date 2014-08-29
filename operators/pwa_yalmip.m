function varargout = pwa_yalmip(varargin)
%PWA_YALMIP Defines a piecewise function using data from MPT
%
%Only intended for internal use in YALMIP
%
% Three classes of PWA functions can be generated
%
% 1) Convexity aware epigraph version with a
%    convex domain and objective modelled as c'x <t.
%    Note, if used in non-concvex setting, a milp
%    model will be generated as a back-up.
%
%    This is the function generated to describe value
%    function of a single mpLP problem.
%
% 2) Overlapping partitions with convex functions
%    in each partition. Requires one binary per region.
%
%    This is the function used to describe the value
%    function generated when not removing overlaps in
%    binary mpLP.
%
% 3) Non-overlapping general pwa function. Requires
%    one binary per region.
%
%    This is the function generated for value functions
%    when remove overlaps is done. It is also used to
%    describe the pwa optimizer.
%
% The input is a cell with structs. Each struct has
% the fields Pn, Pfinal, Bi and Ci (or Fi and Gi)

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
                    if isfield(pwastruct{1},'valuefunction')
                        % Several solutions, pick the one related to the lowest value function
                        valfcn = [];
                        for k = 1:length(jj)
                            valfcn = [valfcn pwastruct{i}.valuefunction.Bi{jj(k)}*x+pwastruct{i}.valuefunction.Ci{jj(k)}];
                        end
                        [iii,jjj] = min(valfcn);
                        val = min(val,pwastruct{i}.Bi{jj(jjj)}(index,:)*x+pwastruct{i}.Ci{jj(jjj)}(index));
                    else
                        % Several solutions, pick the lowest one
                        for k = 1:length(jj)
                            val = min(val,pwastruct{i}.Bi{jj(k)}(index,:)*x+pwastruct{i}.Ci{jj(k)}(index));
                        end
                    end
                end
            end
            if isinf(val)
                val = nan;
            end
            varargout{1} = val;
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
                error('Wrong format on input to PWA (requires Bi or Fi)');
            else
                for i = 1:length(varargin{1})
                    varargin{1}{i}.Bi = varargin{1}{i}.Fi
                    varargin{1}{i}.Ci = varargin{1}{i}.Gi
                end
            end
        end

        if ~isfield(varargin{1}{1},'Pfinal')
            error('Wrong format on input to PWA (requires field Pn)');
        end

        % Create binary variables already here to avoid generating
        % binary variables for the same variable several times
        switch varargin{3}
            case 'convex'
                % No binary needed
                varargin{end+1} = [];
            case 'convexoverlapping'
                % One binary per overlap
                varargin{end+1} = binvar(length(varargin{1}),1);
            case 'general'
                % one binary per region
                varargin{end+1} = binvar(length(varargin{1}{1}.Pn),1);
            otherwise
        end

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

                % Can only generate the first class of PWA functions
                t     = varargin{2};     % The YALMIP variables modelling this pwa
                pwa_struct = varargin{3};% MPT structure
                x     = varargin{4};     % Argument
                flag  = varargin{5};     % Type of PWA function
                d     = varargin{6};     % Binary for nonconvex cases
                index = varargin{7};     % Which row in Bix+Ci

                switch flag

                    case 'convex'

                        [H,K] = double(pwa_struct{1}.Pfinal);
                        S = unique(([reshape([pwa_struct{1}.Bi{:}]',length(x),[])' reshape([pwa_struct{1}.Ci{:}]',[],1)]),'rows');
                       % S = unique(round(S*1e8),'rows')/1e8;
                        
                        costs = S*[x;1];
                        %costs = reshape([pwa_struct{1}.Bi{:}]',length(x),[])'*x+reshape([pwa_struct{1}.Ci{:}]',[],1);
                        F = (H* x <= K) + ((costs<=t):'Epigraph of pwa');

                    case 'convexoverlapping'

                        % Local costs
                        s = sdpvar(length(pwa_struct),1);

                        % Some region, one defined as minimizer
                        F = (sum(d)==1);
                        maxcost = -inf;
                        mincost = inf;
                        for i = 1:length(pwa_struct)
                            % Reduce complexity of max function by
                            % removing equal costs
                            B = reshape([pwa_struct{i}.Bi{:}]',length(x),[])';
                            C = reshape([pwa_struct{i}.Ci{:}]',[],1);
                            [ii,jj,kk] = unique(round(1e8*[B C])/1e8,'rows');
                            cost = reshape([pwa_struct{i}.Bi{jj}]',length(x),[])'*x+reshape([pwa_struct{i}.Ci{jj}]',[],1);
                            [M,m] = derivebounds(cost);
                            maxcost = max(maxcost,max(M));
                            mincost = min(mincost,min(m));
                            [H,K] = double(pwa_struct{i}.Pfinal);
                            [M,m] = derivebounds(H*x - K);
                            F = F + (H*x-K <= (1+M)*(1-d(i)));
                        end
                        for i = 1:length(pwa_struct)

                            B = reshape([pwa_struct{i}.Bi{:}]',length(x),[])';
                            C = reshape([pwa_struct{i}.Ci{:}]',[],1);
                            [ii,jj,kk] = unique(round(1e8*[B C])/1e8,'rows');

                            cost = reshape([pwa_struct{i}.Bi{jj}]',length(x),[])'*x+reshape([pwa_struct{i}.Ci{jj}]',[],1);

                            [M,m] = derivebounds(cost);
                            F = F + (cost <= t + 2*(1+maxcost)*(1-d(i)));
                        end
                        [t_bounds] = yalmip('getbounds',getvariables(t));
                        bounds(t,max([t_bounds(1) mincost]),min([t_bounds(2) 3*maxcost]));


                    case 'general'

                        % In one region
                        F = (sum(d) == 1);

                        % Extract the wanted row
                        for i = 1:length(pwa_struct{1}.Pn)
                            pwa_struct{1}.Bi{i} = pwa_struct{1}.Bi{i}(index,:);
                            pwa_struct{1}.Ci{i} = pwa_struct{1}.Ci{i}(index);
                        end

                        % value in different regions
                        costs = reshape([pwa_struct{1}.Bi{:}]',length(x),[])'*x+reshape([pwa_struct{1}.Ci{:}]',[],1);
                        [M,m] = derivebounds(costs);
                                            
                        % t equals some of the costs
                        % Big-M : |costs-t| < max(costs)-min(t) < M-m
                        F = F + ((m-max(M)).*(1-d)  <= costs-t <= (M-min(m)).*(1-d));
                        for i = 1:length(pwa_struct{1}.Pn)
                            [H,K] = double(pwa_struct{1}.Pn(i));
                            [M,m] = derivebounds(H*x - K);
                            F = F + (H*x - K <= M*(1-d(i)));
                        end

                    otherwise

                        varargout{1} = [];
                        varargout{2} = struct('convexity','milp','monotonicity','none','definiteness','none');
                        varargout{3} = [];
                        return
                end

                varargout{1} = F;
                varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','none','model','graph');
                varargout{3} = x;

            case 'exact'

                % FIX : Should create case for overlapping convex PWAs,
                % used in a nonconvex fashion...

                % Can only generate the first class of PWA functions
                t     = varargin{2};     % The YALMIP variables modelling this pwa
                pwa_struct = varargin{3};% MPT structure
                x     = varargin{4};     % Argument
                flag  = varargin{5};     % Type of PWA function
                d     = varargin{6};     % Binary for nonconvex cases
                index = varargin{7};     % Which row in Bix+Ci
                
                switch flag

                    case {'general','convex'}

                        if isempty(d)
                            d = binvar(length(pwa_struct{1}.Pn),1)
                        end
                        F = (sum(d) == 1);

                        % Extract the wanted row
                        for i = 1:length(pwa_struct{1}.Pn)
                            pwa_struct{1}.Bi{i} = pwa_struct{1}.Bi{i}(index,:);
                            pwa_struct{1}.Ci{i} = pwa_struct{1}.Ci{i}(index);
                        end

                        [p,lb,ub]=bounding_box(pwa_struct{1}.Pn);

                        % Cost in different regions
                        H = reshape([pwa_struct{1}.Bi{:}]',length(x),[])';
                        K = reshape([pwa_struct{1}.Ci{:}]',[],1);
                        costs = H*x+K;
                        [M,m] = derivebounds(costs);
                        M2 = H.*(H>0)*ub + H.*(H<0)*lb+K;
                        m2 = H.*(H>0)*lb + H.*(H<0)*ub+K;
                        M = min([M M2],[],2);
                        m = max([m m2],[],2);

                        % Might be called with a convex function,
                        % but wants the complete MILP description
                        % Example : Someone tries to maximize value
                        % function
                        if length(d)~=length(costs)
                            d = binvar(length(costs),1);
                            F = (sum(d) == 1);
                        end
                        % t equals some of the costs
                        F = F + ((m-max(M)).*(1-d)  <= costs-t <= (M-min(m)).*(1-d));
                        for i = 1:length(pwa_struct{1}.Pn)
                            [H,K] = double(pwa_struct{1}.Pn(i));
                            [M,m] = derivebounds(H*x - K);
                            M2 = H.*(H>0)*ub + H.*(H<0)*lb-K;
                            m2 = H.*(H>0)*lb + H.*(H<0)*ub-K;
                            M = min([M M2],[],2);
                            m = max([m m2],[],2);
                            F = F + (H*x - K <= M*(1-d(i)));
                        end

                    otherwise

                        varargout{1} = [];
                        varargout{2} = struct('convexity','convex','monotonicity','none','definiteness','none','model','integer')
                        varargout{3} = [];
                        return
                end

                varargout{1} = F;
                varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
                varargout{3} = x;


            otherwise
                error('PWA_YALMIP called with CHAR argument?');
        end
    otherwise
        error('Strange type on first argument in SDPVAR/NORM');
end