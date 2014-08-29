function F = sdpvar(X)
% SDPVAR Converts a block variable to standard SDPVAR

% first, check that we have filled everything.
[n,m] = size(X.blocks);

if n==m
    % Perform symmetric extension
    for i = 1:n
        for j = i+1:n
            if isempty(X.blocks{i,j})
                X.blocks{i,j} = X.blocks{j,i}';
            elseif isempty(X.blocks{j,i})
                X.blocks{j,i} = X.blocks{i,j}';
            end
        end
    end
end

for i = 1:n
    for j = 1:m
        ns(i,j) = size(X.blocks{i,j},1);
        ms(i,j) = size(X.blocks{i,j},2);
        % Indicate zero/eye blocks (simplies code below)
        if ns(i,j)==1 & ms(i,j) == 1 & isequal(X.blocks{i,j},0)
            ns(i,j) = 0;
            ms(i,j) = 0;
        end
    end
end

newns = ns;
newms = ms;

% Expand missing elements
for i = 1:n
    for j = 1:m
        if ns(i,j) == 0
            ns(i,j) = max(ns(i,:));
        end
        if ms(i,j) == 0
            ms(i,j) = max(ms(:,j));
        end
    end
end

if any(any(ns==0)) | any(any(ms==0))
    error('I cannot determine the implied size of some of the blocks');
end

F = [];
for i = 1:n
    Y = [];
    for j = 1:m
        if isequal(X.blocks{i,j},0)
            Y = [Y zeros(ns(i,j),ms(i,j))];
        elseif isequal(X.blocks{i,j},1)
            Y = [Y eye(ns(i,j))];
        elseif isempty(X.blocks{i,j}) & isempty(X.blocks{j,i})
            Y = [Y zeros(ns(i,j),ms(i,j))];
        else
            Y = [Y X.blocks{i,j}];
        end
    end
    F = [F;Y];
end