function Z = replace(X,Y,W,expand)
%REPLACE Substitutes variables
%
%Z = REPLACE(X,Y,W)  Replaces any occurence of the SDPVAR object Y
%                    in the SDPVAR object X with the expression W
%
% Example
%  x = sdpvar(1,1);
%  t = sdpvar(1,1);
%  Y = [1+t;1+x+t];
%  Y = replace(Y,x,2) generates Y=[1+t;3+t]

if nargin<4
    expand = 1;
end

if ~isa(X,'sdpvar')
    Z = X;
    return
end
if ~isa(Y,'sdpvar')
    error('Second arguments must be an sdpvar object')
end

if ~is(Y,'linear') 
    error('Second arguments must be linear')
end

if prod(size(W)) == 1
    W = repmat(W,size(Y));
end

if ~isequal(size(Y),size(W))
    if isequal(fliplr(size(Y)),size(W))
        W = W';
    else
        error('Both arguments must have same size')
    end
end

if max(size(Y))>1
 [Y,keptY] = unique(reshape(Y,[],1));
 W = extsubsref(W,keptY);
end

if isa(W,'sdpvar')
    % This is tricky...
    Z = variable_replace(X,Y,W);
    return
end

if ~isnumeric(W)
    error('Third arguments must be numeric')
end

% Replace with NaN   destroys everything, assume it should be cleared
W(isnan(W)) = 0;

y_lmi_variables = Y.lmi_variables;
b = W(:)-Y.basis(:,1);
A = Y.basis(:,2:end);
feas_var = A\b;
if norm(A*feas_var-b)>sqrt(eps)
    error('Inconsistent assignment')
end

x_lmi_variables = X.lmi_variables;
n = X.dim(1);
m = X.dim(2);

[monomtable,variabletype] = yalmip('monomtable');
if all(variabletype(x_lmi_variables)==0) % is(X,'linear')
    Z = X.basis(:,1);
    %v = [];
    v1 = [];
    v2 = [];
    i1 = [];
    i2 = [];
    for i = 1:length(x_lmi_variables)
        j = find(x_lmi_variables(i) == y_lmi_variables);        
        if isempty(j)
          %  v = [v ;recover(x_lmi_variables(i))];
            v1 = [v1 x_lmi_variables(i)];
            i1 = [i1 i];
            %Z = Z + recover(x_lmi_variables(i))*X.basis(:,i+1);
        else
          %  v = [v ;feas_var(j)];
            v2 = [v2 j];
            i2 = [i2 i];
            %Z = Z + feas_var(j)*X.basis(:,i+1);
        end
    end
    v = sparse(i1,ones(length(i1),1),recover(v1),length(x_lmi_variables),1);
    v = v + sparse(i2,ones(length(i2),1),feas_var(v2),length(x_lmi_variables),1);
    Z = Z + X.basis(:,2:end)*v;
else
    base = getbase(Y);base = base(:,2:end);
    [i,j,k] = find(base);
    replaced_vars = getvariables(Y);
    replaced_vars = replaced_vars(i);
    %for i = 1:length(Y)
    %    replaced_vars(i) = getvariables(extsubsref(Y,i));
    %end
    % used_variables = getvariables(X);
    used_variables = x_lmi_variables;
    %  monomtable = yalmip('monomtable');
    local_monom = monomtable(used_variables,replaced_vars);
    W = W(:)';
    gain = zeros(length(used_variables),1);
    for i = 1:length(used_variables)
        % F**N 6.5 0^sparse(0) and 0^0 differ
        gain(i) = prod(W.^full(local_monom(i,:)));
    end

    local_monoms_left = monomtable(used_variables,:);
    local_monoms_left(:,replaced_vars) = 0;
    used_left = find(sum(local_monoms_left,1));
    base = recovermonoms(local_monoms_left(:,used_left),recover(used_left));
    base = base.*gain(:);
    Z = X.basis(:,1);
    Z = Z + X.basis(:,2:end)*base;
end

if expand
    Xvariables = getvariables(Z);
    extvar = yalmip('extvariables');
    Xext = find(ismember(Xvariables,extvar));
    if ~isempty(Xext)
        %We must dig down in extended operators to see if they use the replaced
        %set of variables
        for i = 1:length(Xext)
            extstruct = yalmip('extstruct',Xvariables(Xext(i)));
            anychanged = 0;
            for j = 1:length(extstruct.arg)
                if isa(extstruct.arg{j},'sdpvar')
                    XinY = find(ismembc(getvariables(extstruct.arg{j}),y_lmi_variables));
                    if ~isempty(XinY)
                        anychanged = 1;
                        extstruct.arg{j} = replace(extstruct.arg{j},Y,W);
                    else
                    end
                end
            end
            if anychanged
                if isequal(extstruct.fcn,'pwa_yalmip') | isequal(extstruct.fcn,'pwq_yalmip')
                    % Change data in MPT structure to allow plotting of PWA
                    % and PWQ functions (and improve numerics)
                    % FIXME : Generalize code from sdpvar/plot
                    % extstruct = derive_new_mpt_data(extstruct);
                end
                Zi = yalmip('define',extstruct.fcn,extstruct.arg{:});
                Xvariables(Xext(i)) = getvariables(Zi);
            end
        end
        % And now recover this sucker
        Z = struct(Z);
        Z.lmi_variables = Xvariables;
        % Fucked up order (lmi_variables should be sorted)
        if any(diff(Z.lmi_variables)<0)
            [i,j]=sort(Z.lmi_variables);
            Z.basis = [Z.basis(:,1) Z.basis(:,j+1)];
            Z.lmi_variables = Z.lmi_variables(j);
        end
        Z = sdpvar(Z.dim(1),Z.dim(2),[],Z.lmi_variables,Z.basis);
    end
end

if isa(Z,'sdpvar')
    Z.dim(1) = n;
    Z.dim(2) = m;
    Z.typeflag = X.typeflag;
    Z.extra = X.extra;
    % Reset info about conic terms
    Z.conicinfo = [0 0];
else
    Z = reshape(full(Z),n,m);
end
