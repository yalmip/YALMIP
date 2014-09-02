function [ML,list] = exponents(poly,x)
%EXPONENTS Internal function to extract powers of nonlinear expression

% HOrrible code for noncommuting change. Basically silly change of the
% commuting case, should be written from scratch instead
mt = yalmip('monomtable');
nonCommutingTable = yalmip('nonCommutingTable');
x_lin = getvariables(poly);

x_var = getvariables(x);

ML = 0*mt(x_lin,x_var);
zero = 0;
if any(full(poly.basis(:,1))) %any(ismember(1,poly))
    ML = 0*[zeros(1,length(x));ML];
    zero = 1;
end

list = [];
for i = 1:length(x_lin)
    variable = x_lin(i);
    row = [];
    if nnz(mt(variable,:)) == 0
        % noncommuting term, start with the commuting scaling
        if ~isnan(nonCommutingTable(variable,1))
            scaling = nonCommutingTable(variable,1);
            uses = find(mt(scaling,:));
            for j = 1:length(uses)
                if ismember(uses(j),x_var)
                    row = [row repmat(uses(j),1,mt(scaling,uses(j)))];
                end
            end
        end

        row = [row nonCommutingTable(variable,2:end)];

    else
        uses = find(mt(variable,:));
        for j = 1:length(uses)
            if ismember(uses(j),x_var)
                row = [row repmat(uses(j),1,mt(variable,uses(j)))];
            end
        end
    end
    if isempty(row)
        row = 0;
    end
    row = row(find(row));
    for j = 1:length(row)
        if ismember(row(j),x_var)
            row(j) = find(x_var == row(j));
        else
            row(j) = 0;
        end
    end
    row = row(find(row));
    list(end+1,1:length(row)) = row;
end
ML = zeros(size(list,1),length(x_var));
for i = 1:size(list,1)
    j = 1;
    while j<= size(list,2) & list(i,j)
        ML(i,list(i,j)) =  ML(i,list(i,j))+1;
        j = j + 1;
    end
end

