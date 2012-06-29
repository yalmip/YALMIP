function y = minus(X,Y)
%MINUS (overloaded)

% Author Johan Löfberg
% $Id: minus.m,v 1.2 2006-08-11 11:48:15 joloef Exp $

if isa(X,'sdpvar')
    X = ncvar(struct(X));
elseif isa(Y,'sdpvar')
    Y = ncvar(struct(Y));
end

X_is_ncvar = isa(X,'ncvar');
Y_is_ncvar = isa(Y,'ncvar');

switch 2*X_is_ncvar+Y_is_ncvar
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

        % Speeeeeeed
        if x_isscalar & y_isscalar
             y.basis = -y.basis;
             y.basis(1) = y.basis(1)+X;
             % Reset info about conic terms
             y.conicinfo = [0 0];
             return
        end
         
        if any_scalar | ([n_Y m_Y]==[n_X m_X])
            if y_isscalar
                y.basis = repmat(y.basis,n_X*m_X,1);
                y.dim(1) = n_X;
                y.dim(2) = m_X;
            end
            y.basis = -y.basis;
            if nnz(X)~=0
                y.basis(:,1) = y.basis(:,1)+X(:);
            end
        else
            error('Matrix dimensions must agree.');
        end
        % Reset info about conic terms
        y.conicinfo = [0 0];

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

        % Silly hack
        % Taking X-scalar(0) takes unnecessary time
        % and is used in most definitions of LMIs
        if (y_isscalar & (Y==0))
            return
        end
                
        % Speeeeeeed
        if x_isscalar & y_isscalar
             y.basis(1) = y.basis(1)-Y;
             % Reset info about conic terms
             y.conicinfo = [0 0];
             return
         end

        if any_scalar | ([n_Y m_Y]==[n_X m_X])
            if x_isscalar
                y.basis = repmat(y.basis,n_Y*m_Y,1);
                y.dim(1) = n_Y;
                y.dim(2) = m_Y;
            end
            y.basis(:,1) = y.basis(:,1)-Y(:);
        else
            error('Matrix dimensions must agree.');
        end

        % Update information about conic terms
        % This information is used in DUALIZE to 
        % speed up some checks, and to facilitate some
        % advanced dualization features. It also
        % speeds up checking for symmetry in some other code
        % Ugly, but the best way at the moment
        % For a description of this field, check SDPVAR code
       % if (y.conicinfo(1)~=0) & isequal(Y,Y') & (y.conicinfo(2) ~= 2)
       %     y.conicinfo(2) = max(1,y.conicinfo(2));
       % else
            y.conicinfo = [0 0];
       % end
        

    case 3

        %	if (X.typeflag~=0) | (Y.typeflag~=0)
        %		error('Relational objects cannot be manipulated')
        %	end

        n_X = X.dim(1);
        m_X = X.dim(2);
        n_Y = Y.dim(1);
        m_Y = Y.dim(2);
        x_isscalar = (n_X*m_X==1);
        y_isscalar = (n_Y*m_Y==1);
        any_scalar = x_isscalar | y_isscalar;

        if ~any_scalar
            if (~((n_X==n_Y) & (m_X==m_Y)))
                error('Matrix dimensions must agree.')
            end
        end

        all_lmi_variables = uniquestripped([X.lmi_variables Y.lmi_variables]);
        y = X;
        X.basis = []; % Returns memory?
        y.lmi_variables = all_lmi_variables;

        in_X_logical = ismembc(all_lmi_variables,X.lmi_variables);
        in_Y_logical = ismembc(all_lmi_variables,Y.lmi_variables);
        in_X = find(in_X_logical);
        in_Y = find(in_Y_logical);

        if isequal(X.lmi_variables,Y.lmi_variables) & n_Y==n_X & m_Y==m_X
            y.basis = y.basis - Y.basis;
            % Super special case f(scalar)-f(scalar)
            if length(X.lmi_variables)==1
                if all(y.basis(:,2)==0)
                    y = full(y.basis(1));
                else
                    y.conicinfo = [0 0];
                end
                return
            end
        else
            if 1
                [ix,jx,sx] = find(y.basis);y.basis = [];
                [iy,jy,sy] = find(Y.basis);Y.basis = [];
                mapX = [1 1+in_X];
                mapY = [1 1+in_Y];
                basis_X = sparse(ix,mapX(jx),sx,n_X*m_X,1+length(all_lmi_variables));ix=[];jx=[];sx=[];
                basis_Y = sparse(iy,mapY(jy),sy,n_Y*m_Y,1+length(all_lmi_variables));iy=[];jy=[];sy=[];
            else
                % MATLAB sparse fails on this for huge problems at a certain size
                basis_X = spalloc(n_X*m_X,1+length(all_lmi_variables),nnz(X.basis));
                basis_Y = spalloc(n_Y*m_Y,1+length(all_lmi_variables),nnz(Y.basis));
                basis_X(:,[1 1+in_X])=y.basis;y.basis = [];
                basis_Y(:,[1 1+in_Y])=Y.basis;Y.basis = [];
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
                basis_Y = repmat(basis_Y,n_X*m_X,1);
            end

            % OK, solution is...
            y.basis = basis_X - basis_Y;
        end
        
        % Only clean if there are variables used in both
        %if ~all(xor(in_X_logical,in_Y_logical))
        % Reset info about conic terms
                 
     %    if (y.conicinfo(1)~=0) & ishermitian(Y) & isempty(intersect(X.lmi_variables,Y.lmi_variables))
     %       y.conicinfo = [y.conicinfo(2);
     %   else
            y.conicinfo = [0 0];
     %   end
         
         
         y = clean(y);
        %else
        %end

    otherwise
end

% Update info on KYP objects
if X_is_ncvar & Y_is_ncvar & X.typeflag==9  & Y.typeflag==9
    error('Substraction of KYP objects currently not supported')
end
if Y_is_ncvar & Y.typeflag==9
    y.extra.M = -Y.extra.M+X;
    y.extra.negated = ~Y.extra.negated;
    return
end
if X_is_ncvar & X.typeflag==9
    y.extra.M = y.extra.M-Y;
    return
end



