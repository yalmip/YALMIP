function [cuts,maxViolation] = separateGramTriangles(gram,x,nvars)
% Cut a fractional point with inequalities certified by 3x3 PSD minors.
% A pivot j represents the two equal edges (i,j) and (j,k).
rows = [];
columns = [];
entries = [];
count = 0;
maxViolation = 0;
violations = [];
for block = 1:numel(gram)
    pair = gram{block}.pair;
    equalValue = gram{block}.equalValue;
    n = size(pair,1);
    valid = pair > 0;
    value = zeros(n);
    value(valid) = x(pair(valid));
    if equalValue == 1
        constant = 1;
        coefficients = [1 -1 -1];
    else
        constant = 0;
        coefficients = [-1 1 1];
    end
    for pivot = 1:n
        side = bsxfun(@plus,value(:,pivot),value(pivot,:));
        if equalValue == 1
            violation = side-value-1;
        else
            violation = value-side;
        end
        usable = bsxfun(@and,valid(:,pivot),valid(pivot,:)) & valid;
        violation(~usable) = 0;
        violation = triu(violation,1);
        maxViolation = max(maxViolation,max(violation(:)));
        [i,k] = find(violation > 1e-6);
        m = numel(i);
        if m == 0
            continue
        end
        opposite = pair(sub2ind([n n],i,k));
        first = pair(i,pivot);
        second = pair(pivot,k)';
        cutRows = (count+(1:m))';
        rows = [rows;repmat(cutRows,4,1)]; %#ok<AGROW>
        columns = [columns;ones(m,1);1+opposite;1+first;1+second]; %#ok<AGROW>
        entries = [entries;constant*ones(m,1); ...
            coefficients(1)*ones(m,1);coefficients(2)*ones(m,1); ...
            coefficients(3)*ones(m,1)]; %#ok<AGROW>
        violations = [violations;violation(sub2ind([n n],i,k))]; %#ok<AGROW>
        count = count+m;
    end
end
cuts = sparse(rows,columns,entries,count,nvars+1);
if count > 5000
    [~,order] = sort(violations,'descend');
    cuts = cuts(order(1:5000),:);
end
end
