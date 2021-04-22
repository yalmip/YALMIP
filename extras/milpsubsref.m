function varargout = milpsubsref(varargin)
%MILPSUBSREF

switch class(varargin{1})
    case 'double'
        varargin{1} = double(varargin{1});
        for i = 1:length(varargin{2}.subs)
            varargin{2}.subs{i} = round(double(varargin{2}.subs{i}));
        end
        varargout{1} = subsref(varargin{:});

    case {'sdpvar','ndsdpvar'}
        switch length(varargin{2}.subs)
            case 1
                index = varargin{2}.subs{1};
                if length(index) > 1
                    y = [];
                    for i = 1:length(index)
                        varargin{2}.subs{1} = index(i);
                        y = [y;yalmip('define',mfilename,varargin{:})];
                    end
                    % Figure out dims of variable
                    X = randn(size(varargin{1}));
                    varargin{2}.subs{1} = ones(size(index));
                    X = subsref(X,varargin{2});
                    y = reshape(y,size(X,1),size(X,2));
                    varargout{1} = y;
                else
                    varargout{1} = yalmip('define',mfilename,varargin{:});
                end
            case 2
                index1 = varargin{2}.subs{1};
                if isa(index1,'char')
                    index1 = 1:size(varargin{1},1);
                end
                index2 = varargin{2}.subs{2};
                if isa(index2,'char')
                    index2 = 1:size(varargin{1},2);
                end
                y = [];
                if length(index1)*length(index2) == 1
                    varargin{2}.subs{1} = index1;
                    varargin{2}.subs{2} = index2;
                    y = yalmip('define',mfilename,varargin{:});
                else
                    for i = 1:length(index1)
                        temp = [];
                        for j = 1:length(index2)
                            varargin{2}.subs{1} = index1(i);
                            varargin{2}.subs{2} = index2(j);
                            temp = [temp yalmip('define',mfilename,varargin{:})];
                        end
                        y = [y;temp];
                    end
                end
                % Figure out dims of variable
                X = randn(size(varargin{1}));
                varargin{2}.subs{1} = ones(size(index1));
                varargin{2}.subs{2} = ones(size(index2));
                X = subsref(X,varargin{2});
                y = reshape(y,size(X,1),size(X,2));
                varargout{1} = y;
                
            case 3
                index1 = varargin{2}.subs{1};
                if isa(index1,'char')
                    index1 = 1:size(varargin{1},1);
                end
                index2 = varargin{2}.subs{2};
                if isa(index2,'char')
                    index2 = 1:size(varargin{1},2);
                end
                index3 = varargin{2}.subs{3};
                if isa(index3,'char')
                    index3 = 1:size(varargin{1},3);
                end
                y = [];
                if length(index1)*length(index2)*length(index3) == 1
                    varargin{2}.subs{1} = index1;
                    varargin{2}.subs{2} = index2;
                    varargin{2}.subs{3} = index3;
                    y = yalmip('define',mfilename,varargin{:});
                else
                    temp = [];
                    for i = 1:length(index1)                        
                        for j = 1:length(index2)
                            for k = 1:length(index3)
                                varargin{2}.subs{1} = index1(i);
                                varargin{2}.subs{2} = index2(j);
                                varargin{2}.subs{3} = index3(k);
                                temp = [temp yalmip('define',mfilename,varargin{:})];
                            end
                        end                        
                    end                    
                end
                % Figure out dims of variable
                X = randn(size(varargin{1}));
                varargin{2}.subs{1} = ones(size(index1));
                varargin{2}.subs{2} = ones(size(index2));
                varargin{2}.subs{3} = ones(size(index3));
                X = subsref(X,varargin{2});
                y = reshape(y,size(X,1),size(X,2),size(X,3));
                varargout{1} = y;                
            otherwise
                error('Only 1D and 2D variable subsref implemented');
        end

    case 'char'

        X = varargin{3};
        Y = varargin{2};
        if length(varargin{4}.subs) == 1
            X = X(:);
            i = varargin{4}.subs{1};
            M = length(X);
            m = 1;
            F = (integer(i)); % just to be sure
            d = binvar(length(X),1);
            [Mx,mx]=derivebounds(X);
            for j = m:M
                di = d(j);
                %                        F = F + (mx*(1-di) <= Y-X(j) <= Mx*(1-di));
                F = F + (-(max(Mx)-min(mx))*(1-di) <= Y-X(j) <= (max(Mx)-min(mx))*(1-di));
                F = F + (-(1+M-m)*(1-di) <= i-j <= (1+M-m)*(1-di));
            end
            F = F + (sum(d)==1);
        elseif length(varargin{4}.subs) == 2
            i1 = varargin{4}.subs{1};
            i2 = varargin{4}.subs{2};
            M1 = size(X,1);
            M2 = size(X,2);
            m1 = 1;
            m2 = 1;
            if isa(i1,'sdpvar')
                F = (integer(i1)); % just to be sure
            end
            if isa(i2,'sdpvar')
                F = (integer(i2)); % just to be sure
            end
            d = binvar(size(X,1),size(X,2),'full');
            [Mx,mx]=derivebounds(X);
            for i = m1:M1
                for j = m2:M2
                    di = d(i,j);
                    F = F + (-(max(Mx)-min(mx))*(1-di) <= Y-X(i,j) <= (max(Mx)-min(mx))*(1-di));
                    if isa(i1,'sdpvar')
                        F = F + (-(1+M1-m1)*(1-di) <= i1-i <= (1+M1-m1)*(1-di));
                    elseif i1~=i
                        F = [F, sum(d(i,:))==0];
                    end
                    if isa(i2,'sdpvar')
                        F = F + (-(1+M2-m2)*(1-di) <= i2-j <= (1+M2-m2)*(1-di));
                    elseif i2~=j
                        F = [F, sum(d(:,j))==0];
                    end
                end
            end
            F = F + (sum(sum(d))==1);
        elseif length(varargin{4}.subs) == 3
            i1 = varargin{4}.subs{1};
            i2 = varargin{4}.subs{2};
            i3 = varargin{4}.subs{3};
            M1 = size(X,1);
            M2 = size(X,2);
            M3 = size(X,3);
            m1 = 1;
            m2 = 1;
            m3 = 1;
            if isa(i1,'sdpvar')
                F = (integer(i1)); % just to be sure
            end
            if isa(i2,'sdpvar')
                F = (integer(i2)); % just to be sure
            end
            if isa(i3,'sdpvar')
                F = (integer(i3)); % just to be sure
            end
            d = binvar(size(X,1),size(X,2),size(X,3),'full');
            [Mx,mx]=derivebounds(X(:));            
            for i = m1:M1
                for j = m2:M2
                    for k = m3:M3
                        di = d(i,j,k);
                        F = F + (-(max(Mx)-min(mx))*(1-di) <= Y-X(i,j,k) <= (max(Mx)-min(mx))*(1-di));
                        if isa(i1,'sdpvar')
                            F = F + (-(1+M1-m1)*(1-di) <= i1-i <= (1+M1-m1)*(1-di));
                        elseif i1~=i
                            F = [F, sum(d(i,:,:))==0];
                        end
                        if isa(i2,'sdpvar')
                            F = F + (-(1+M2-m2)*(1-di) <= i2-j <= (1+M2-m2)*(1-di));
                        elseif i2~=j
                            F = [F, sum(d(:,j,:))==0];
                        end
                        if isa(i3,'sdpvar')
                            F = F + (-(1+M3-m3)*(1-di) <= i3-k <= (1+M3-m3)*(1-di));
                        elseif i3~=k
                            F = [F, sum(d(:,:,k))==0];
                        end
                    end
                end
            end
            F = F + (sum(sum(sum(d)))==1);                                                            
        end
        varargout{1} = F;
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = [X(:);i(:)];

    otherwise
        error('Strange type on first argument in SDPVAR/MILPSUBSREF');
end
