function Z = ncvar_replace(X,Y,W)

% Call recurisvely for vector X
if length(X) > 1
    dimX = size(X);
    Z = [];
    X = X(:);
    for i = 1:length(X)
        Z = [Z;ncvar_replace(X(i),Y,W)];
    end
    Z = reshape(Z,dimX);
    return
end

% ...and recursively for vector Y
if length(Y) > 1    
    Z = X;
    for i = 1:length(Y)
        Z = ncvar_replace(Z,Y(i),W(i));
    end
    return
end

nonCommutingTable         = yalmip('nonCommutingTable');
[monomtable,variabletype] = yalmip('monomtable');
if size(monomtable,1)>size(nonCommutingTable,1)
    nonCommutingTable((1+size(nonCommutingTable,1)):(size(monomtable,1)),1) = (1+size(nonCommutingTable,1)):(size(monomtable,1));
end

% Cast commutative variables as nc temporarily by adding them to the table
commuting = find(~any(nonCommutingTable,2));
nonCommutingTable(commuting,1) = commuting;

% Expand non-commuting nonlinears
for i = 1:size(nonCommutingTable)
    if variabletype(i) & ~isnan(nonCommutingTable(i,1))
        [aux,monoms,powers] = find(monomtable(i,:));
        r = [];
        for j = 1:length(monoms)
            r = [r repmat(monoms(j),1,powers(j))];
        end
        nonCommutingTable(i,1:length(r))=r;
    end
end

base = getbase(X);
vars = getvariables(X);
yvar = getvariables(Y);
Z = X;
newMonoms = [1];
for i = vars
    monTerms = nonCommutingTable(i,:);
    monTerms = monTerms(monTerms ~= 0);
    monTerms = monTerms(~isnan(monTerms));
    if any(ismember(yvar,monTerms))
        newMonomial = 1;
        for j = 1:length(monTerms)
            if monTerms(j)==yvar
                newMonomial = newMonomial*W;
            else
                newMonomial = newMonomial*recover(monTerms(j));
            end
        end
        newMonoms = [newMonoms;newMonomial];
    else
        newMonoms = [newMonoms;recover(i)];
    end
end
Z = reshape(base*newMonoms,size(X));