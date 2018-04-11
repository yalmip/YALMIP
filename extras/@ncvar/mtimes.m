function Z = mtimes(X,Y)
%MTIMES (overloaded)

% Brute-force implementation of multiplication of noncommuting variables
% (with possible commuting variables involved)

% Check classes
X_is_spdvar = isa(X,'sdpvar');
Y_is_spdvar = isa(Y,'sdpvar');

X_is_ncvar = isa(X,'ncvar');
Y_is_ncvar = isa(Y,'ncvar');

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

x_variables = getvariables(X);Xbase = getbase(X);
y_variables = getvariables(Y);Ybase = getbase(Y);
temp_monom_table = [];
temp_nc_table = [];
temp_c_table = [];
new_base = [];

for i = 0:length(x_variables)
    if i>0
        x_monom = nonCommutingTable(x_variables(i),:);
    else
        x_monom = nan;
    end
    x_base = Xbase(:,i+1);
    for j = 0:length(y_variables)
        if j>0
            y_monom = nonCommutingTable(y_variables(j),:);
        else
            y_monom = nan;
        end
        y_base = Ybase(:,j+1);
        xy_base = reshape(x_base,size(X))*reshape(y_base,size(Y));

        if (i == 0) & (j== 0)
            new_base = xy_base(:);
        elseif nnz(xy_base)>0
            xy_monom = [x_monom(2:end) y_monom(2:end)];
            xy_monom = xy_monom(find(xy_monom));
            temp_nc_table(end+1,1:length(xy_monom)) = xy_monom;
            temp_c_table(end+1,1:2) = [x_monom(1) y_monom(1)];
            new_base = [new_base xy_base(:)];
        end
    end
end

% It could have happended that new commuting monomials where generated
% during the multiplication. Check and create these
for i = 1:size(temp_c_table,1)
    aux = spalloc(1,size(monomtable,2),2);
    if ~isnan(temp_c_table(i,1))
        aux = monomtable(temp_c_table(i,1),:) + aux;
    end
    if ~isnan(temp_c_table(i,2))
        aux = monomtable(temp_c_table(i,2),:)  + aux;
    end
    if nnz(aux)>0
        candidates = findrows(monomtable,aux);
        if ~isempty(candidates)
            temp_c_table(i,1) = candidates;
        else
            monomtable = [monomtable;aux];
            nonCommutingTable(end+1,1) = nan;
            temp_c_table(i,1) = size(monomtable,1);
            switch sum(aux)
                case 1
                    variabletype(end+1) = 0;
                case 2
                    if nnz(aux) == 1
                        variabletype(end+1) = 2;
                    else
                        variabletype(end+1) = 1;
                    end
                otherwise
                    variabletype(end+1) = 3;
            end
        end
    end
end

temp_nc_table = [temp_c_table(:,1) temp_nc_table];
% Okay, now we have the monomials. Now we have to match them to
% possible earlier monomials
if size(nonCommutingTable,2) < size(temp_nc_table,2)
    nonCommutingTable(1,size(temp_nc_table,2)) = 0;
elseif size(temp_nc_table,2) < size(nonCommutingTable,2)
    temp_nc_table(1,size(nonCommutingTable,2)) = 0;
end
for i = 1:size(temp_nc_table,1)
    candidates = findrows_nan(nonCommutingTable,temp_nc_table(i,:));
    if isempty(candidates)
        nonCommutingTable = [nonCommutingTable;temp_nc_table(i,:)];
        monomtable(end+1,end+1) = 0;
        involved =  temp_nc_table(i,1+find(temp_nc_table(i,2:end)));
        switch length(involved)
            case 1
                if isnan(temp_nc_table(i,1))
                    variabletype(end+1) = 0;
                else
                    if variabletype(temp_nc_table(i,1)) == 0
                        variabletype(end+1) = 1;
                    else
                        variabletype(end+1) = 3;
                    end
                end
            case 2
                if isnan(temp_nc_table(i,1))
                    if involved(1) == involved(2)
                        variabletype(end+1) = 2;
                    else
                        variabletype(end+1) = 1;
                    end
                else
                    variabletype(end+1) = 3;
                end
            otherwise
                variabletype(end+1) = 3;
        end
        lmivariables(i) = size(nonCommutingTable,1);
    else
        lmivariables(i) = candidates;

    end
end

% Create an output variable quickly
if X_is_ncvar
    Z = X;
else
    Z = Y;
end
Z.basis = new_base;
Z.lmi_variables = lmivariables;

% Fucked up order (lmi_variables should be sorted and unique)
if any(diff(Z.lmi_variables)<0)
    [i,j]=sort(Z.lmi_variables);
    Z.basis = [Z.basis(:,1) Z.basis(:,j+1)];
    Z.lmi_variables = Z.lmi_variables(j);
end
[un_Z_vars2] = uniquestripped(Z.lmi_variables);
if length(un_Z_vars2) < length(Z.lmi_variables)
    [un_Z_vars,hh,jj] = unique(Z.lmi_variables);
    if length(Z.lmi_variables) ~=length(un_Z_vars)
        Z.basis = Z.basis*sparse([1 1+jj(:)'],[1 1+(1:length(jj))],ones(1,1+length(jj)))';
        Z.lmi_variables = un_Z_vars;
    end
end
Z.dim = size(xy_base);
Z = clean(Z);
if size(monomtable,2) < size(monomtable,1)
    monomtable(size(monomtable,1),size(monomtable,1)) = 0;
end
yalmip('nonCommutingTable',nonCommutingTable);
yalmip('setmonomtable',monomtable,variabletype);

function c = findrows_nan(a,b)
a(isnan(a)) = 0;
b(isnan(b)) = 0;
c=findrows(a,b);