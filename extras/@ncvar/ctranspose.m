function X=ctranspose(X)
%CTRANSPOSE (overloaded)

if isa(X,'blkvar')
    X = sdpvar(X);
end

n = X.dim(1);
m = X.dim(2);
ind = reshape(reshape(1:n*m,n,m)',n*m,1);
if isreal(X.basis)
    X.basis = X.basis(ind,:);    
else
   X.basis = conj(X.basis(ind,:));
end
X.dim(1) = m;
X.dim(2) = n;
% Reset info about conic terms
X.conicinfo = [0 0];

% Flip noncommuting terms
% Get all the tables, and expand them so that they correspond to the same
% number of variables globally (nonCommutingTable is not up to date after a
% new commuting variables has been defined, to save flops)
nonCommutingTable         = yalmip('nonCommutingTable');
[monomtable,variabletype] = yalmip('monomtable');
if size(monomtable,1)>size(nonCommutingTable,1)
    nonCommutingTable((1+size(nonCommutingTable,1)):(size(monomtable,1)),1) = (1+size(nonCommutingTable,1)):(size(monomtable,1));
end
% Cast commutative variables as nc temporarily by adding them to the table
commuting = find(~any(nonCommutingTable,2));
nonCommutingTable(commuting,1) = commuting;

for i = 1:length(X.lmi_variables)
    if nonCommutingTable(X.lmi_variables(i),2)
        monoms = nonCommutingTable(X.lmi_variables(i),2:end);      
        monoms(find(monoms)) = fliplr( monoms(find(monoms)));
        monoms = [nonCommutingTable(X.lmi_variables(i),1) monoms];
        old = findrows(nonCommutingTable,monoms);
        old = findrows(nonCommutingTable(:,2:end),monoms(2:end));
        if ~isempty(old)
            old = old(find((monoms(1) == nonCommutingTable(old,1)) | isnan(monoms(1)) & isnan(nonCommutingTable(old,1))));
        end
        if isempty(old)
            % Create a new monomial
            monomtable(end+1,end+1) = 0;
            variabletype(end+1) = variabletype(X.lmi_variables(i));
            nonCommutingTable = [nonCommutingTable;monoms];
            X.lmi_variables(i) = size(nonCommutingTable,1);
        else
            X.lmi_variables(i) = old;
        end
    end
end
% Fucked up order (lmi_variables should be sorted and unique)
if any(diff(X.lmi_variables)<0)
    [i,j]=sort(X.lmi_variables);
    X.basis = [X.basis(:,1) X.basis(:,j+1)];
    X.lmi_variables = X.lmi_variables(j);
end
[un_Z_vars2] = uniquestripped(X.lmi_variables);
if length(un_Z_vars2) < length(X.lmi_variables)
    [un_Z_vars,hh,jj] = unique(X.lmi_variables);
    if length(X.lmi_variables) ~=length(un_X_vars)
        X.basis = Z.basis*sparse([1 1+jj],[1 1+(1:length(jj))],ones(1,1+length(jj)))';
        X.lmi_variables = un_Z_vars;
    end
end
if size(monomtable,2) < size(monomtable,1)
    monomtable(size(monomtable,1),size(monomtable,1)) = 0;
end
yalmip('nonCommutingTable',nonCommutingTable);
yalmip('setmonomtable',monomtable,variabletype);
