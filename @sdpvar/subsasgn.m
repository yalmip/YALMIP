function y = subsasgn(X,I,Y)
%SUBASGN (overloaded)

try
    if strcmp('()',I.type)
        X_is_spdvar = isa(X,'sdpvar') |  isa(X,'ndsdpvar');
        Y_is_spdvar = isa(Y,'sdpvar') |  isa(Y,'ndsdpvar');
        if islogical(I.subs{1})
            I.subs{1} = double(find(I.subs{1}));
        end
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
                    X0 = subsasgn(full(X),I,full(reshape(Y.basis(:,1),n_y,m_y)));
                    dim = size(X0);
                    y.basis = reshape(X0,prod(dim),1);
                    X = full(X)*0;
                    for i = 1:length(y_lmi_variables)
                        X0 = subsasgn(X,I,full(reshape(Y.basis(:,i+1),n_y,m_y)));
                        y.basis(:,i+1) = reshape(X0,prod(dim),1);
                    end
                    y.dim = dim;
                    % Reset info about conic terms
                    y.conicinfo = [0 0];
                    y.basis = sparse(y.basis);
                    if length(dim)>2
                        y = ndsdpvar(y);
                    end
                    y = flush(y);
                catch
                    error(lasterr)
                end
            case 2
                if ~isempty(Y)     
                    if isa(Y,'uint8') || isa(Y,'uint16') || isa(Y,'uint32') || isa(Y,'uint64')
                        Y = sparse(double(Y));
                    elseif isnumeric(Y)
                        Y = sparse(Y);
                    else
                        Y = sparse(double(Y));
                    end
                end
                y = X;
                
                % Special code for speed
                % elements in vector replaced with constants
                if min(X.dim(1),X.dim(2))==1 & (length(I.subs)==1)
                     y = X;
                     if isempty(Y)
                         y.basis(I.subs{1},:) = [];
                         if X.dim(1) == 1
                             y.dim(2) = y.dim(2) - length(unique(I.subs{1}));
                         else
                             y.dim(1) = y.dim(1) - length(unique(I.subs{1}));
                         end
                     else
                         y.basis(I.subs{1},1) = Y;
                         y.basis(I.subs{1},2:end) = 0;                         
                     end
                     if prod(y.dim)~=size(y.basis,1)
                         % Ah bugger, the dimension of the object was changed)
                         aux = X.basis(:,1);
                         aux = reshape(aux,X.dim);
                         aux(I.subs{1})=Y;
                         y.dim = size(aux);
                     end
                     y = clean(y);
                     % Reset info about conic terms
                     if isa(y,'sdpvar')
                         y.conicinfo = [0 0];
                         y = flush(y);
                     end
                     return;
                end
                    
                
                x_lmi_variables = X.lmi_variables;
                lmi_variables = [];
                           
                n = y.dim(1);
                m = y.dim(2);
                subX = sparse(subsasgn(full(reshape(X.basis(:,1),n,m)),I,Y));
                y.basis = subX(:);
                if isa(I.subs{1},'char')
                    I.subs{1} = 1:n;
                end
                if length(I.subs)>1
                    if isa(I.subs{2},'char')
                        I.subs{2} = 1:m;
                    end
                end
                if length(I.subs)>1
                    if length(I.subs{1})==1 & length(I.subs{2})~=1
                        I.subs{1} = repmat(I.subs{1},size(I.subs{2},1),size(I.subs{2},2));
                    elseif length(I.subs{2})==1 & length(I.subs{1})~=1
                        I.subs{2} = repmat(I.subs{2},size(I.subs{1},1),size(I.subs{1},2));
                    end
                end
                
                if length(I.subs)>1                  
                    ii = kron(I.subs{1}(:),ones(length(I.subs{2}),1));
                    jj = kron(ones(length(I.subs{1}),1),I.subs{2}(:));
                    LinearIndex = sub2ind([n m],ii,jj);
                else
                    LinearIndex = I.subs{1};
                end
                
                if isempty(Y)
                    X.basis = X.basis(:,2:end);
                    X.basis(LinearIndex,:) = [];
                    y.basis = [y.basis(:,1) X.basis];
                else
                    X.basis(LinearIndex,2:end)=sparse(0);                
                    y.basis = [y.basis(:,1) X.basis(:,2:end)];
                end
                         
                y.dim(1) = size(subX,1);
                y.dim(2) = size(subX,2);
                y = clean(y);
                if isa(y,'sdpvar')
                    % Reset info about conic terms
                    y.conicinfo = [0 0];
                    y = flush(y);
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
                
                if (mx==1) & (my == 1) & isempty(setdiff(y_lmi_variables,x_lmi_variables)) & (max(I.subs{1}) < nx) & length(I.subs)==1 & length(unique(I.subs{1}))==length(I.subs{1}) ;
                    % Fast specialized code for Didier
                     y = specialcode(X,Y,I);
                     return
                end                
               
                subX = subsasgn(reshape(X.basis(:,1),nx,mx),I,reshape(Y.basis(:,1),ny,my));
                [newnx, newmx] = size(subX);
                                              
                j = 1;
                
                yz = reshape(1:ny*my,ny,my);
                subX2 = subsasgn(reshape(zeros(nx*mx,1),nx,mx),I,yz);
                subX2 = subX2(:);
                [ix,jx,sx] = find(subX2);
                yz = 0*reshape(Y.basis(:,1),ny,my);                               
                lmi_variables = zeros(1,length(x_lmi_variables));
                
                A = reshape(1:nx*mx,nx,mx);
                B = reshape(1:newnx*newmx,newnx,newmx);
                
                rm = B(1:nx,1:mx);rm = rm(:);
                [iix,jjx,ssx] = find(X.basis(:,2:end));
                z.basis = [subX(:) sparse(rm(iix),jjx,ssx,newnx*newmx,size(X.basis,2)-1)];
                z.basis(ix,2:end) = 0;
                                               
                keep = find(any(z.basis(:,2:end),1));
                z.basis = z.basis(:,[1 1+keep]);
                lmi_variables2 = x_lmi_variables(keep);

                z.lmi_variables = lmi_variables2;
                lmi_variables = lmi_variables2;
                
                all_lmi_variables = union(lmi_variables,y_lmi_variables);
                in_z = ismembcYALMIP(all_lmi_variables,lmi_variables);
                in_y = ismembcYALMIP(all_lmi_variables,y_lmi_variables);
                z_ind = 2;
                y_ind = 2;
                basis = spalloc(size(z.basis,1),1+length(all_lmi_variables),0);
                basis(:,1) = z.basis(:,1);
                nz = size(subX,1);
                mz = size(subX,2);
                template = full(0*reshape(X.basis(:,1),nx,mx));
                in_yin_z = 2*in_y + in_z;
                if all(in_yin_z<3)
                    case1 = find(in_yin_z==1);
                    if ~isempty(case1)
                        basis(:,case1+1) = z.basis(:,2:1+length(case1));                        
                        in_yin_z(case1) = 0;
                    end
                end
                    % Let's identify the indices at which we need to intervene
                    case1 = find(in_yin_z==1);
                    checkI = union(find(in_yin_z >= 2), setdiff(case1,case1+1));
                    if size(checkI,1) > size(checkI,2)
                        % Sometimes the in_yin_z vector is vertical (not
                        % always though), so we make sure that checkI, on
                        % which we'll iterate, is horizontal.
                        checkI = checkI';
                    end
                    for i = checkI
                        switch in_yin_z(i)
                            case 1
                                % We look for the end of the block of ones
                                % starting at i
                                iend = i + find([in_yin_z(i+1:end) 0] ~= 1, 1, 'first')-1;
                                basis(:,i+1:iend+1) = z.basis(:,z_ind:z_ind+(iend-i));z_ind = z_ind+(iend-i)+1;
                        case 2                          
                            temp = sparse(subsasgn(template,I,full(reshape(Y.basis(:,y_ind),ny,my))));
                            basis(:,i+1) = temp(:);
                            y_ind = y_ind+1;
                        case 3
                            Z1 = z.basis(:,z_ind);
                            Z4 = Y.basis(:,y_ind);
                            Z3 = reshape(Z4,ny,my);
                            Z2 = sparse(subsasgn(0*reshape(full(X.basis(:,1)),nx,mx),I,Z3));
                            temp = reshape(Z1,nz,mz)+Z2;                            
                            basis(:,i+1) = temp(:);
                            z_ind = z_ind+1;
                            y_ind = y_ind+1;
                        otherwise
                    end 
                end;
                z.dim(1) = nz;
                z.dim(2) = mz;
                z.basis = basis;
                z.lmi_variables = all_lmi_variables(:)';
                y = z;	                
                % Reset info about conic terms
                y.conicinfo = [0 0];                 
                y = flush(y);
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


