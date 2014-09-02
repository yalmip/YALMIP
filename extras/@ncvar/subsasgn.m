function y = subsasgn(X,I,Y)
%SUBASGN (overloaded)

try
    if strcmp('()',I.type)
        X_is_spdvar = isa(X,'sdpvar');
        Y_is_spdvar = isa(Y,'sdpvar');
        if any(I.subs{1} <=0)
            error('Index into matrix is negative or zero.');
        end
        switch 2*X_is_spdvar+Y_is_spdvar
            case 1 
                % This code does not work properly
                % Only work if b is undefined!!?!!
                % generally ugly code...
                y = Y;
                [n_y,m_y] = size(Y);
                y_lmi_variables = y.lmi_variables;
                try
                    X0 = sparse(subsasgn(full(X),I,full(reshape(Y.basis(:,1),n_y,m_y))));
                    [n_x,m_x] = size(X0);
                    y.basis = reshape(X0,n_x*m_x,1);
                    X = full(X)*0;
                    for i = 1:length(y_lmi_variables)
                        X0 = full(sparse(subsasgn(X,I,full(reshape(Y.basis(:,i+1),n_y,m_y)))));
                        y.basis(:,i+1) = reshape(X0,n_x*m_x,1);
                    end
                    y.dim(1) = n_x;
                    y.dim(2) = m_x;
                    % Reset info about conic terms
                    y.conicinfo = [0 0];
                catch
                    error(lasterr)
                end
            case 2
                if ~isempty(Y)
                  Y = sparse(Y);
                end
                y = X;
                
                % Special code for speed
                % elements in vector replaced with constants
                if min(X.dim(1),X.dim(2))==1 & (length(I.subs)==1)
                     y = X;
                     y.basis(I.subs{1},1) = Y;
                     y.basis(I.subs{1},2:end) = 0;                     
                     y = clean(y);
                     % Reset info about conic terms
                     if isa(y,'sdpvar')
                         y.conicinfo = [0 0];
                     end
                     return;
                end
                    
                
                x_lmi_variables = X.lmi_variables;
                lmi_variables = [];
                
                % y.basis = [];
                n = y.dim(1);
                m = y.dim(2);
                subX = sparse(subsasgn(full(reshape(X.basis(:,1),n,m)),I,Y));
                y.basis = subX(:);
                
                j = 1;
                Z = 0*Y;
                for i = 1:length(x_lmi_variables)
                    subX = sparse(subsasgn(full(reshape(X.basis(:,i+1),n,m)),I,Z));
                    if (norm(subX,inf)>0)
                        y.basis(:,j+1) = subX(:);
                        lmi_variables = [lmi_variables x_lmi_variables(i)];
                        j = j+1;
                    end
                end  
                y.dim(1) = size(subX,1);
                y.dim(2) = size(subX,2);
                if isempty(lmi_variables) % Convert back to double!!
                    y=full(reshape(y.basis(:,1),y.dim(1),y.dim(2)));
                    return
                else %Nope, still a sdpvar
                    y.lmi_variables = lmi_variables;
                     % Reset info about conic terms
                    y.conicinfo = [0 0];
                end
                
            case 3
                z = X;
                
                x_lmi_variables = X.lmi_variables;
                y_lmi_variables = Y.lmi_variables;
                
                                
                % In a first run, we fix the constant term and null terms in the X basis
                lmi_variables = [];
                nx = X.dim(1);
                mx = X.dim(2);
                ny = Y.dim(1);
                my = Y.dim(2);
                
                if (mx==1) & (my == 1) & isempty(setdiff(y_lmi_variables,x_lmi_variables)) & (max(I.subs{1}) < nx);
                    % Fast specialized code for Didier
                     y = specialcode(X,Y,I);
                     return
                end
                %      subX = sparse(subsasgn(full(reshape(X.basis(:,1),nx,mx)),I,reshape(Y.basis(:,1),ny,my)));
                subX = subsasgn(reshape(X.basis(:,1),nx,mx),I,reshape(Y.basis(:,1),ny,my));
                
                z.basis = subX(:);
                j = 1;
                yz = 0*reshape(Y.basis(:,1),ny,my);
                for i = 1:length(x_lmi_variables)
                    %        subX = sparse(subsasgn(full(reshape(X.basis(:,i+1),nx,mx)),I,yz));
                    subX = subsasgn(reshape(X.basis(:,i+1),nx,mx),I,yz);
                    
                    if (norm(subX,inf)>0)
                        z.basis(:,j+1) = subX(:);
                        lmi_variables = [lmi_variables x_lmi_variables(i)];
                        j = j+1;
                    end
                end
                z.lmi_variables=lmi_variables;
                all_lmi_variables = union(lmi_variables,y_lmi_variables);
                in_z = ismembc(all_lmi_variables,lmi_variables);
                in_y = ismembc(all_lmi_variables,y_lmi_variables);
                z_ind = 2;
                y_ind = 2;
                basis=z.basis(:,1);
                nz = size(subX,1);
                mz = size(subX,2);
                for i = 1:length(all_lmi_variables)
                    switch 2*in_y(i)+in_z(i)
                        case 1
                            basis(:,i+1) = z.basis(:,z_ind);z_ind = z_ind+1;
                        case 2
                            temp = sparse(subsasgn(full(0*reshape(X.basis(:,1),nx,mx)),I,full(reshape(Y.basis(:,y_ind),ny,my))));
                            basis(:,i+1) = temp(:);
                            y_ind = y_ind+1;
                        case 3
                            Z1 = z.basis(:,z_ind);
                            Z4 = Y.basis(:,y_ind);
                            Z3 = reshape(Z4,ny,my);
                            Z2 = sparse(subsasgn(0*reshape(full(X.basis(:,1)),nx,mx),I,Z3));
                            temp = reshape(Z1,nz,mz)+Z2;                            
%                            temp = reshape(z.basis(:,z_ind),nz,mz)+sparse(subsasgn(0*reshape(full(X.basis(:,1)),nx,mx),I,reshape(Y.basis(:,y_ind),ny,my)));
                            basis(:,i+1) = temp(:);
                            z_ind = z_ind+1;
                            y_ind = y_ind+1;
                        otherwise
                    end 
                end;
                z.dim(1) = nz;
                z.dim(2) = mz;
                z.basis = basis;
                z.lmi_variables = all_lmi_variables;
                y = z;	                
                % Reset info about conic terms
                y.conicinfo = [0 0];                 
            otherwise
        end
    else
        error('Reference type not supported');
    end
    
catch
    error(lasterr)
end


function y = specialcode(X,Y,I)

y = X;
X_basis = X.basis;
Y_basis = Y.basis;
ind = I.subs{1};ind = ind(:);
yvar_in_xvar = zeros(length(Y.lmi_variables),1);
for i = 1:length(Y.lmi_variables);
    yvar_in_xvar(i) = find(X.lmi_variables==Y.lmi_variables(i));
end
y.basis(ind,:) = 0;
mapper = [1 1+yvar_in_xvar(:)'];mapper = mapper(:);
[i,j,k] = find(y.basis);
[ib,jb,kb] = find(Y_basis);
i = [i(:);ind(ib(:))];
j = [j(:);mapper(jb(:))];
k = [k(:);kb(:)];
y.basis = sparse(i,j,k,size(y.basis,1),size(y.basis,2));
y = clean(y);


