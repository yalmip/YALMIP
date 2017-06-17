function Z = addfactors(Z,X,Y)
if isnumeric(X) || isa(X,'logical')
    if length(Y.midfactors)==0
        return
    end
    %   Y = refactor(Y);
    Z.midfactors = Y.midfactors;
    Z.leftfactors = Y.leftfactors;
    Z.rightfactors = Y.rightfactors;
    Z.midfactors{end+1} = X;
    if length(X)>1
        Z.leftfactors{end+1} = speye(size(X,1));
        Z.rightfactors{end+1} = speye(size(X,2));
    else
        Z.leftfactors{end+1} = 1;
        Z.rightfactors{end+1} = 1;
    end
elseif isnumeric(Y) || isa(Y,'logical')
    if length(X.midfactors)==0
        return
    end
    %    X = refactor(X);
    Z.midfactors = X.midfactors;
    Z.leftfactors = X.leftfactors;
    Z.rightfactors = X.rightfactors;
    Z.midfactors{end+1} = Y;
    if length(Y)>1
        Z.leftfactors{end+1} = speye(size(Y,1));
        Z.rightfactors{end+1} = speye(size(Y,2));
    else
        Z.leftfactors{end+1} = 1;
        Z.rightfactors{end+1} = 1;
    end
else
    if length(X.midfactors)==0 || length(Y.midfactors)==0
        Z.midfactors = [];
        Z.leftfactors = [];
        Z.rightfactors = [];
        return
    end
    %    X = refactor(X);
    %    Y = refactor(Y);
    
    if prod(X.dim)>0 &&  prod(Y.dim)==1
        for i = 1:length(Y.midfactors)
            Y.leftfactors{i} = repmat(Y.leftfactors{i},X.dim(1),1);
            Y.rightfactors{i} = repmat(Y.rightfactors{i},1,X.dim(2));
        end
    end
      if prod(Y.dim)>0 &&  prod(X.dim)==1
        for i = 1:length(X.midfactors)
            X.leftfactors{i} = repmat(X.leftfactors{i},Y.dim(1),1);
            X.rightfactors{i} = repmat(X.rightfactors{i},1,Y.dim(2));
        end
    end
    
    Z.leftfactors = {X.leftfactors{:},Y.leftfactors{:}};
    Z.midfactors = {X.midfactors{:},Y.midfactors{:}};
    Z.rightfactors = {X.rightfactors{:},Y.rightfactors{:}};
end
n = Z.dim(1);
m = Z.dim(2);

for i = 1:length(Z.midfactors)
   % isdouble(i) = isa(Z.midfactors{i},'double');
    if size(Z.leftfactors{i},1)~=n
        if prod(size(Z.midfactors{i}))==1
            Z.leftfactors{i} = ones(n,1);
        else
            error('Inconsistent factors. Please report bug');
        end
    end
    if size(Z.rightfactors{i},2)~=m
        if prod(size(Z.midfactors{i}))==1
            Z.rightfactors{i} = ones(1,m);
        else
            error('Inconsistent factors. Please report bug');
        end
    end
end
Z = cleandoublefactors(Z);
try
Z = flushmidfactors(Z);
catch
    1
end

function Z = refactor(Z)
if isempty(Z.midfactors)
    Z.leftfactors{1} = 1;
    Z.midfactors{1} = Z;
    Z.rightfactors{1} = 1;
end
