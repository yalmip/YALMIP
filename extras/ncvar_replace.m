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

base = getbase(X);
vars = getvariables(X);
yvar = getvariables(Y);
Z = X;
newMonoms = [1];
for i = vars
    monTerms = nonCommutingTable(i,:);
    
    % Possibilities
    % [v(x) NC]
    % [NAN NC]
    if isa(Y,'sdpvar')       
        if isnan(monTerms(1))
            % Purely non-commuting so no replacement can occur
            newMonoms = [newMonoms;recover(i)];
        else
            % call recursively with sdpvar replace
            if any(isnan(nonCommutingTable(monTerms(find(monTerms)))))
                [ct,vt] = coefficients(recovernc(i),Y);
            else
                [ct,vt] = coefficients(recover(i),Y);
            end
            newMonomial = ct;
            if degree(vt)>=1
                newMonomial = newMonomial*W^degree(vt);
            end
            newMonoms = [newMonoms;newMonomial];
        end
    else
        % Replacing a non-commuting with something else
        newMonomial = 1;
        for j = 2:max(find(monTerms))
            if monTerms(j)==yvar
                newMonomial = newMonomial*W;
            else
                newMonomial = newMonomial*recover(monTerms(j));
            end
        end
        newMonoms = [newMonoms;newMonomial];
    end                 
end
Z = reshape(base*newMonoms,size(X));
