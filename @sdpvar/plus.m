function y = plus(X,Y)
%PLUS (overloaded)

% Cannot use isa here since blkvar is marked as sdpvar
X_class = class(X);
Y_class = class(Y);
X_is_spdvar = strcmp(X_class,'sdpvar');
Y_is_spdvar = strcmp(Y_class,'sdpvar');
    
% Convert block objects
if ~X_is_spdvar
    if isa(X,'blkvar')
        X = sdpvar(X);
        X_is_spdvar = isa(X,'sdpvar');
    elseif isa(X,'intval')
        X_is_spdvar = 0;
        Y.basis = intval(Y.basis);
    elseif isa(X,'uint8') || isa(X,'uint16') || isa(X,'uint32') || isa(X,'uint64')
        X = double(X);  
    elseif ~isnumeric(X)
        error(['Cannot add SDPVAR object and ' upper(class(X)) ' object']);
    end
end

if ~Y_is_spdvar
    if isa(Y,'blkvar')
        Y = sdpvar(Y);
        Y_is_spdvar = isa(Y,'sdpvar');;
    elseif isa(Y,'intval')
        Y_is_spdvar = 0;
        X.basis = intval(X.basis);
    elseif isa(Y,'uint8') || isa(Y,'uint16') || isa(Y,'uint32') || isa(Y,'uint64')
        Y = double(Y);          
    elseif ~isnumeric(Y)
        error(['Cannot add SDPVAR object and ' upper(class(Y)) ' object']);
    end
end

if isnumeric(X)
    if any(isnan(X))
        disp('You have NaNs in model (<a href="yalmip.github.io/naninmodel">learn to debug</a>)')
        error('Adding NaN to an SDPVAR makes no sense.');
    end
end
if isnumeric(Y)
     if any(isnan(Y))
        disp('You have NaNs in model (<a href="yalmip.github.io/naninmodel">learn to debug</a>)')
        error('Adding NaN to an SDPVAR makes no sense.');
     end
end

switch 2*X_is_spdvar+Y_is_spdvar
    case 1
        if isempty(X)
            try
                y = full(X - reshape(Y.basis(:,1),Y.dim(1),Y.dim(2)));
            catch
                error(lasterr);
            end
            return
        end

        y = Y;
        n_Y = Y.dim(1);
        m_Y = Y.dim(2);
        [n_X,m_X] = size(X);
        x_isscalar = (n_X*m_X==1);
        y_isscalar = (n_Y*m_Y==1);
        any_scalar = x_isscalar | y_isscalar;
        
        if x_isscalar && y_isscalar            
            tmp = y.basis(1)+X;
            if (isequal(class(tmp),'gem') || isequal(class(tmp),'sgem')) && ~isequal(class(y.basis), class(tmp))
                y.basis = gemify(y.basis);
            end
            y.basis(1) = tmp;
            % Reset info about conic terms
            y.conicinfo = [0 0];
            y.extra.opname='';
            y.extra.createTime = definecreationtime;            
            return
         end
         
        if any_scalar || all([n_Y m_Y]==[n_X m_X])
            if y_isscalar
                y.basis = repmat(y.basis,n_X*m_X,1);
                y.dim(1) = n_X;
                y.dim(2) = m_X;
            end
            tmp = y.basis(:,1)+X(:);
            if (isequal(class(tmp),'gem') || isequal(class(tmp),'sgem')) && ~isequal(class(y.basis), class(tmp))
                y.basis = gemify(y.basis);
            end
            y.basis(:,1) = tmp;
        else
            error('Matrix dimensions must agree.');
        end
        % Reset info about conic terms
        y.conicinfo = [0 0];
        y.extra.opname='';        
    case 2

        if isempty(Y)
            try
                y = full(reshape(X.basis(:,1),X.dim(1),X.dim(2))-Y);
            catch
                error(lasterr);
            end
            return
        end

        y = X;
        n_X = X.dim(1);
        m_X = X.dim(2);
        [n_Y,m_Y] = size(Y);
        x_isscalar = (n_X*m_X==1);
        y_isscalar = (n_Y*m_Y==1);
        any_scalar = x_isscalar | y_isscalar;
        
         % Special special case...
         if x_isscalar && y_isscalar
             tmp = y.basis(1)+Y;
             if (isequal(class(tmp),'gem') || isequal(class(tmp),'sgem')) && ~isequal(class(y.basis), class(tmp))
                 y.basis = gemify(y.basis);
             end
             y.basis(1) = tmp;
             % Reset info about conic terms
             y.conicinfo = [0 0];
             y.extra.opname='';
             y.extra.createTime = definecreationtime;             
             return
         end
         
        if any_scalar || all(([n_Y m_Y]==[n_X m_X]))
            if x_isscalar
                y.basis = repmat(y.basis,n_Y*m_Y,1);
                y.dim(1) = n_Y;
                y.dim(2) = m_Y;
            end
            tmp = y.basis(:,1)+Y(:);
            if (isequal(class(tmp),'gem') || isequal(class(tmp),'sgem')) && ~isequal(class(y.basis), class(tmp))
                y.basis = gemify(y.basis);
            end
            y.basis(:,1) = tmp;
        else
            error('Matrix dimensions must agree.');
        end
        % Reset info about conic terms
        y.conicinfo = [0 0];
        y.extra.opname='';
        y.extra.createTime = definecreationtime;        
             
    case 3

        n_X = X.dim(1);
        m_X = X.dim(2);
        n_Y = Y.dim(1);
        m_Y = Y.dim(2);
        x_isscalar = (n_X*m_X==1);
        y_isscalar = (n_Y*m_Y==1);
        any_scalar = x_isscalar | y_isscalar;

        if (~((n_X==n_Y) && (m_X==m_Y))) && ~any_scalar
            error('Matrix dimensions must agree.')
        end

        if isequal(X.lmi_variables,Y.lmi_variables)
             all_lmi_variables = X.lmi_variables;
             in_X_logical = ones(1,length(all_lmi_variables));
             in_Y_logical = ones(1,length(all_lmi_variables));
        else
            if isempty(Y.lmi_variables)
                all_lmi_variables = [X.lmi_variables];
                in_X_logical = [ones(1,length(X.lmi_variables))];
                in_Y_logical = [zeros(1,length(X.lmi_variables))];
            elseif isempty(X.lmi_variables)
                all_lmi_variables = [Y.lmi_variables];
                in_Y_logical = [ones(1,length(Y.lmi_variables))];
                in_X_logical = [zeros(1,length(Y.lmi_variables))];
            elseif X.lmi_variables(end) < Y.lmi_variables(1)
                all_lmi_variables = [X.lmi_variables Y.lmi_variables];
                in_X_logical = [ones(1,length(X.lmi_variables)) zeros(1,length(Y.lmi_variables))];
                in_Y_logical = [zeros(1,length(X.lmi_variables)) ones(1,length(Y.lmi_variables))];
            elseif X.lmi_variables(1) > Y.lmi_variables(end)
                all_lmi_variables = [Y.lmi_variables X.lmi_variables];
                in_X_logical = [zeros(1,length(Y.lmi_variables)) ones(1,length(X.lmi_variables))];
                in_Y_logical = [ones(1,length(Y.lmi_variables)) zeros(1,length(X.lmi_variables))];
            else
                all_lmi_variables = uniquestripped([X.lmi_variables Y.lmi_variables]);
                in_X_logical = ismembcYALMIP(all_lmi_variables,X.lmi_variables);
                in_Y_logical = ismembcYALMIP(all_lmi_variables,Y.lmi_variables);
            end
        end
        y = X;
        y.lmi_variables = all_lmi_variables;

        % ismembc faster (buggy?)
        in_X = find(in_X_logical);
        in_Y = find(in_Y_logical);

        if isequal(X.lmi_variables,Y.lmi_variables) && n_Y==n_X && m_Y==m_X
            y.basis = y.basis + Y.basis;
             if length(X.lmi_variables)==1
                 if all(y.basis(:,2)==0)
                     y = reshape(full(y.basis(:,1)),n_X,m_X);
                 else
                     % Reset info about conic terms
                     y.conicinfo = [0 0];
                     y.extra.opname='';                     
                 end
                return
            end
        else
            if 1
                if  ~isempty(X.lmi_variables) && ~isempty(Y.lmi_variables) && max(X.lmi_variables) < min(Y.lmi_variables) && n_Y==n_X && m_Y==m_X
                    % special case to speed up Lundback's code massivly
                    % Addition of expressions sharing no variables, with
                    % variables in specific sorted order
                    basis = [y.basis  Y.basis];
                    basis(:,1) = basis(:,1) + Y.basis(:,1);
                    basis(:,size(y.basis,2)+1) = [];
                    y.basis = basis;
                    y.conicinfo = [0 0];
                    y.extra.opname='';                    
                    return
                else
                   % [ix,jx,sx] = find(y.basis);y.basis = [];
                   % [iy,jy,sy] = find(Y.basis);%Y.basis = [];
                    mapX = [1 1+in_X];
                    mapY = [1 1+in_Y];
                   % basis_X = sparse(ix,mapX(jx),sx,n_X*m_X,1+length(all_lmi_variables));ix=[];jx=[];sx=[];
                   % basis_Y = sparse(iy,mapY(jy),sy,n_Y*m_Y,1+length(all_lmi_variables));iy=[];jy=[];sy=[];
                    basis_X = X.basis*(sparse(1:length(mapX),mapX,1,size(X.basis,2),length(all_lmi_variables)+1));
                    basis_Y = Y.basis*(sparse(1:length(mapY),mapY,1,size(Y.basis,2),length(all_lmi_variables)+1));
                end
            else
                % MATLAB sparse fails on this for huge problems at a certain size
                basis_X = spalloc(n_X*m_X,1+length(all_lmi_variables),nnz(X.basis));
                basis_Y = spalloc(n_Y*m_Y,1+length(all_lmi_variables),nnz(Y.basis));
                basis_X(:,[1 1+in_X])=y.basis;y.basis = [];
                basis_Y(:,[1 1+in_Y])=Y.basis;%Y.basis = [];
            end

            % Fix addition of matrix+scalar
            if n_X*m_X<n_Y*m_Y
                y.dim(1) = n_Y;
                y.dim(2) = m_Y;
                basis_X = repmat(basis_X,n_Y*m_Y,1);
            end
            if n_Y*m_Y<n_X*m_X
                y.dim(1) = n_X;
                y.dim(2) = m_X;              
                y.basis = basis_X;basis_X = [];
                try                        
                    y.basis = bsxfun(@plus,y.basis,basis_Y);basis_Y = [];
                catch
                    basis_Y = repmat(basis_Y,n_X*m_X,1);
                    y.basis = y.basis + basis_Y;basis_Y = [];
                end
            else
                % OK, solution is...
                y.basis = basis_X;basis_X = [];
                y.basis = y.basis+basis_Y;basis_Y = [];               
            end
        end
        % Reset info about conic terms
        y.conicinfo = [0 0];
        y.extra.opname=''; 
        y.extra.createTime = definecreationtime;        
        if nnz(in_X_logical & in_Y_logical)>0
            y = clean(y);
        end

    otherwise
end


