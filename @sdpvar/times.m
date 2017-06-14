function y = times(X,Y)
%TIMES (overloaded)

% Check dimensions
[n,m]=size(X);
if ~((prod(size(X))==1) || (prod(size(Y))==1))
    if ~((n==size(Y,1) && (m ==size(Y,2))))
        error('Matrix dimensions must agree.')
    end
end;

% Convert block objects
if isa(X,'blkvar')
    X = sdpvar(X);
end

if isa(Y,'blkvar')
    Y = sdpvar(Y);
end

if isnumeric(X)
    if any(isnan(X))
        error('Multiplying NaN with an SDPVAR makes no sense.');
    end
end
if isnumeric(Y)
    if any(isnan(Y))
        error('Multiplying NaN with an SDPVAR makes no sense.');
    end
end

if isempty(X)
    YY = full(reshape(Y.basis(:,1),Y.dim(1),Y.dim(2)));
    y = X.*YY;
    return
elseif isempty(Y)
    XX = full(reshape(X.basis(:,1),X.dim(1),X.dim(2)));
    y = XX.*Y;
    return
end

if (isa(X,'sdpvar') && isa(Y,'sdpvar'))
    X = flush(X);
    Y = flush(Y);
    if (X.typeflag==5) && (Y.typeflag==5)
        error('Product of norms not allowed');
    end
    
  

    try
        
        y = check_for_special_case(Y,X);
        if ~isempty(y)
            return
        end
        
        % Check for the case x.*y where x and y are unit variables
        [mt,variable_type,hashedMT,hash] = yalmip('monomtable');
        if length(X.lmi_variables)==numel(X) 
            if length(Y.lmi_variables) == numel(Y) 
                if numel(X)==numel(Y)
                    
                    % This looks promising. write as (x0+X)*(y0+Y)
                    X0 = reshape(X.basis(:,1),X.dim);
                    Y0 = reshape(Y.basis(:,1),Y.dim);
                    Xsave = X;
                    Ysave = Y;
                    X.basis(:,1)=0;
                    Y.basis(:,1)=0;
                                          
                    if nnz(X.basis)==numel(X)
                        if nnz(Y.basis)==numel(Y)
                            D = [spalloc(numel(Y),1,0) speye(numel(Y))];
                            if isequal(X.basis,D)
                                if isequal(Y.basis,D)
                                    % Pew.
                                    Z = X;                                    
                                    generated_monoms = mt(X.lmi_variables,:) +  mt(Y.lmi_variables,:);
                                    generated_hash = generated_monoms*hash;
                                    keep = zeros(1,numel(X));
                                    
                                    if all(generated_hash) &&  all(diff(sort([generated_hash;hashedMT])))
                                        Z.lmi_variables =  size(mt,1)+(1:numel(X));
                                        keep = keep + 1;
                                    else
                                        for i = 1:numel(X)
                                            if generated_hash(i)
                                                before = find(abs(hashedMT-generated_hash(i))<eps);
                                                if isempty(before)
                                                    %  mt = [mt;generated_monoms(i,:)];
                                                    keep(i) = 1;
                                                    Z.lmi_variables(i) = size(mt,1)+nnz(keep);
                                                else
                                                    Z.lmi_variables(i) = before;
                                                end
                                            else
                                                Z.lmi_variables(i) = 0;
                                            end
                                        end
                                    end
                                    if any(keep)
                                        keep = find(keep);
                                        mt = [mt;generated_monoms(keep,:)];
                                        yalmip('setmonomtable',mt,[],[hashedMT;generated_hash(keep)],hash);
                                    end
                                                                                                             
                                    if any(diff(Z.lmi_variables)<0)
                                        [i,j]=sort(Z.lmi_variables);
                                        Z.lmi_variables = Z.lmi_variables(j);
                                        Z.basis(:,2:end) = Z.basis(:,j+1);
                                    end
                                    
                                    if Z.lmi_variables(1)==0
                                        i = find(Z.lmi_variables == 0);
                                        Z.basis(:,1) = sum(Z.basis(:,1+i),2);
                                        Z.basis(:,1+i)=[];
                                    end
                                    
                                    Z.conicinfo = [0 0];
                                    Z.extra.opname='';
                                    Z = Z + X0.*Y0 + X0.*Y + X.*Y0;
                                    Z = flush(Z);
                                    y = clean(Z);
                                    return
                                end
                            end
                        end
                    end
                end
                X = Xsave;
                Y = Ysave;
            end
        end
        
        
        x_isscalar =  (X.dim(1)*X.dim(2)==1);
        y_isscalar =  (Y.dim(1)*Y.dim(2)==1);

        all_lmi_variables = uniquestripped([X.lmi_variables Y.lmi_variables]);
        Z = X;Z.dim(1) = 1;Z.dim(2) = 1;Z.lmi_variables = all_lmi_variables;Z.basis = [];

        % Awkward code due to bug in ML6.5
        Xbase = reshape(X.basis(:,1),X.dim(1),X.dim(2));
        Ybase = reshape(Y.basis(:,1),Y.dim(1),Y.dim(2));
        if x_isscalar
            Xbase = sparse(full(Xbase));
        end
        if y_isscalar
            Ybase = sparse(full(Ybase));
        end


        index_Y = zeros(length(all_lmi_variables),1);
        index_X = zeros(length(all_lmi_variables),1);
        for j = 1:length(all_lmi_variables)
            indexy = find(all_lmi_variables(j)==Y.lmi_variables);
            indexx = find(all_lmi_variables(j)==X.lmi_variables);
            if ~isempty(indexy)
                index_Y(j) = indexy;
            end
            if ~isempty(indexx)
                index_X(j) = indexx;
            end
        end

        ny = Y.dim(1);
        my = Y.dim(2);
        nx = X.dim(1);
        mx = X.dim(2);

        % Linear terms
        base = Xbase.*Ybase;
        Z.basis = base(:);
        x_base_not_zero = nnz(Xbase)>0;
        y_base_not_zero = nnz(Ybase)>0;
        for i = 1:length(all_lmi_variables)
            base = 0;
            if index_Y(i) && x_base_not_zero
                base = Xbase.*getbasematrixwithoutcheck(Y,index_Y(i));
            end
            if index_X(i) && y_base_not_zero
                base = base + getbasematrixwithoutcheck(X,index_X(i)).*Ybase;
            end
            Z.basis(:,i+1) = base(:);
        end

        % Nonlinear terms
        i = i+1;
        ix=1;

        new_mt = [];
        %mt = yalmip('monomtable');
        nvar = length(all_lmi_variables);
        local_mt = mt(all_lmi_variables,:);
        theyvars = find(index_Y);
        thexvars = find(index_X);

        hash = randn(size(mt,2),1);
        mt_hash = mt*hash;

        for ix = thexvars(:)'

            %            if mx==1
            Xibase = X.basis(:,1+index_X(ix));
            %           else
            %              Xibase = reshape(X.basis(:,1+index_X(ix)),nx,mx);
            %         end
            mt_x = local_mt(ix,:);

            y_basis = Y.basis(:,1+index_Y(theyvars));
            x_basis = repmat(Xibase,1,length(theyvars(:)'));
            if y_isscalar && ~x_isscalar
                y_basis = repmat(y_basis,nx*mx,1);
            elseif  x_isscalar && ~y_isscalar
                x_basis = repmat(x_basis,ny*my,1);
            end
            allBase = x_basis.*y_basis;
            jjj = 1;
            usedatall = find(any(allBase,1));
            % for iy = theyvars(:)'
            for iy = theyvars(usedatall(:))'
                %  ff=Y.basis(:,1+index_Y(iy));
                %  Yibase = reshape(ff,ny,my);
                %  prodbase = Xibase.*Yibase;
                prodbase = allBase(:,usedatall(jjj));jjj = jjj+1;
                %  prodbase = reshape(prodbase,ny,my);
                if (norm(prodbase,inf)>1e-12)
                    mt_y = local_mt(iy,:);
                    % Idiot-hash the lists
                    new_hash = (mt_x+mt_y)*hash;
                    if abs(new_hash)<eps%if nnz(mt_x+mt_y)==0
                        Z.basis(:,1) = Z.basis(:,1) + prodbase(:);
                    else
                        before = find(abs(mt_hash-(mt_x+mt_y)*hash)<eps);
                        if isempty(before)
                            mt = [mt;mt_x+mt_y];
                            mt_hash = [mt_hash;(mt_x+mt_y)*hash];
                            Z.lmi_variables = [Z.lmi_variables size(mt,1)];
                        else
                            Z.lmi_variables = [Z.lmi_variables before];
                        end
                        Z.basis(:,i+1)  = prodbase(:);i = i+1;
                    end
                end
            end
        end

        % Fucked up order
        if any(diff(Z.lmi_variables)<0)
            [i,j]=sort(Z.lmi_variables);
            Z.lmi_variables = Z.lmi_variables(j);
            Z.basis(:,2:end) = Z.basis(:,j+1);
        end

        % FIX : Speed up
        if length(Z.lmi_variables) ~=length(unique(Z.lmi_variables))
            un_Z_vars = unique(Z.lmi_variables);
            newZbase = Z.basis(:,1);
            for i = 1:length(un_Z_vars)
                newZbase = [newZbase sum(Z.basis(:,find(un_Z_vars(i)==Z.lmi_variables)+1),2)];
            end
            Z.basis = newZbase;
            Z.lmi_variables = un_Z_vars;
        end

        yalmip('setmonomtable',mt);

        if ~(x_isscalar || y_isscalar)
            Z.dim(1) = X.dim(1);
            Z.dim(2) = Y.dim(2);
        else
            Z.dim(1) = max(X.dim(1),Y.dim(1));
            Z.dim(2) = max(X.dim(2),Y.dim(2));
        end
    catch
        error(lasterr)
    end
    % Reset info about conic terms
    Z.conicinfo = [0 0];
    Z.extra.opname='';
    Z = flush(Z);
    y = clean(Z);
    return
end

if isa(X,'sdpvar')
    Z = X;
    X = Y;
    Y = Z;
end

y = Y;
if prod(Y.dim)==1
    y.basis = X(:)*(y.basis);
    y.dim = size(X);
else 
    y.basis = [(Y.basis.')*diag(sparse(X(:)))].';
end

% Reset info about conic terms
y.conicinfo = [0 0];
Z.extra.opname='';
y = flush(y);
y = clean(y);



function y = check_for_special_case(Y,X);
y = [];

if (min(size(X))>1) || (min(size(Y))>1)
    return
end

if ~all(size(Y)==size(X))
    return
end

entropies = zeros(length(Y),1);
if is(X,'linear')
    argst = yalmip('getarguments',Y);
    if length(argst)~=length(X)
        return
    end
    if length(argst) == 1
        args{1} = argst;
    else
        args = argst;
    end
    for i = 1:length(args)
        if isempty(args{i})
            return
        end
        if isequal(args{i}.fcn,'log')
            S(1).subs={i};
            S(1).type='()';
            Z = subsref(X,S);
            if isequal(Z.basis,args{i}.arg{1}.basis)
                if isequal(Z.lmi_variables,args{i}.arg{1}.lmi_variables)
                    entropies(i) = 1;
                end
            end
        end
    end
end
if all(entropies)
    y = -ventropy(X);
end


