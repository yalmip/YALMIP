function gram = initGramTriangleSeparator(p,minSize)
% Recognize binary Gram blocks with a uniform two-value off-diagonal pattern.
% Fixed off-diagonals are allowed, but triangles involving them are skipped.
if nargin<2, minSize=25; end
gram = {};
top = startofSDPCone(p.K);
for cone = 1:numel(p.K.s)
    n = p.K.s(cone);
    % A root cut round has not repaid its work on smaller blocks. Limit the
    % cubic separator work on very large blocks until a cheaper oracle exists.
    if n < minSize || n > 80
        top = top+n*n;
        continue
    end
    block = p.F_struc(top:top+n*n-1,:);
    diagonal = block(1:n+1:n*n,:);
    d = full(diagonal(1,1));
    if nnz(diagonal(:,2:end)) || any(abs(full(diagonal(:,1))-d) > 1e-9)
        top = top+n*n;
        continue
    end
    pair = zeros(n);
    base = NaN;
    coefficient = NaN;
    valid = true;
    exactPattern = all(full(diagonal(:,1)) == d);
    for j = 2:n
        for i = 1:j-1
            row = block(i+(j-1)*n,:);
            mirror = block(j+(i-1)*n,:);
            if norm(row-mirror,inf) > 1e-8
                valid = false;
                break
            end
            if nnz(row-mirror)
                exactPattern = false;
            end
            variable = find(row(2:end));
            if isempty(variable)
                continue
            end
            if numel(variable) ~= 1 || ~p.isbinary(variable)
                valid = false;
                break
            end
            if isnan(base)
                base = full(row(1));
                coefficient = full(row(1+variable));
            elseif abs(full(row(1))-base) > 1e-8 || ...
                    abs(full(row(1+variable))-coefficient) > 1e-8
                valid = false;
                break
            end
            if full(row(1)) ~= base || ...
                    full(row(1+variable)) ~= coefficient
                exactPattern = false;
            end
            pair(i,j) = variable;
            pair(j,i) = variable;
        end
        if ~valid
            break
        end
    end
    if valid && ~isnan(base)
        [~,equalValue] = min(abs([base,base+coefficient]-d));
        equalValue = equalValue-1;
        equalEntry = base+equalValue*coefficient;
        otherEntry = base+(1-equalValue)*coefficient;
        witness = [d equalEntry otherEntry; ...
                   equalEntry d equalEntry; ...
                   otherEntry equalEntry d];
        % A triangle with two equal-value edges and one opposite edge must
        % violate this principal minor. The margin covers SDPA rounding.
        if all(isfinite(witness(:))) && ...
                min(eig(witness)) < -1e-5*(1+norm(witness,'fro'))
            gram{end+1} = struct('pair',pair,'equalValue',equalValue, ...
                'saturated',equalEntry == d,'diagonal',d, ...
                'otherEntry',otherEntry,'exactPattern',exactPattern); %#ok<AGROW>
        end
    end
    top = top+n*n;
end
end
