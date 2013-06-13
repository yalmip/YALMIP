function [c,v] = ncvar_coefficients(p,x)

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

base = getbase(p);
vars = getvariables(p);
xvar = getvariables(x);
c = base(1);
v = [1];
for i = vars
    monTerms = nonCommutingTable(i,:);
    monTerms = monTerms(monTerms ~= 0);
    monTerms = monTerms(~isnan(monTerms));
    if any(ismember(xvar,monTerms))
        newMonomial = 1;
        newBase = 1;
        for j = 1:length(monTerms)
            if any(monTerms(j)==xvar)
                newMonomial = newMonomial*recover(monTerms(j));
            else
                newBase = newBase*recover(monTerms(j));
            end
        end
        c = [c;newBase];
        v = [v;newMonomial];       
    else
        c = c + recover(i);
    end
end
if isequal(c(1),0)
    v = v(2:end);
    c = c(2:end);
end