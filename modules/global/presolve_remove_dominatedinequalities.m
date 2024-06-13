function p = presolve_remove_dominatedinequalities(p)
% Simple removal of b2 >= a'*x if b1 >= a'*x with b1 < b2 exist
if p.K.l > 0 && p.feasible
    % Lazy, so I start by sorting them
    bA = p.F_struc(1+p.K.f:p.K.f+p.K.l,:);
    [~,sortidx] = sort(bA(:,1),'ascend');
    bA = bA(sortidx,:);
    hash = bA*[0;randn(length(p.c),1)];
    [~,idx] = unique(hash,'stable');
    if length(idx) < size(bA,1)
        keepthese = sortidx(idx);
        removethese = p.K.f + setdiff(1:p.K.l,keepthese);
        p.F_struc(removethese,:) = [];
        n = length(removethese);
        p.K.l = p.K.l - n;
        if p.options.verbose > 1 && length(n)>0
            disp(['* Removed ' num2str(n) ' dominated inequalities']);
        end
    end
end