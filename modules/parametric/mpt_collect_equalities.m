function Matrices = mpt_collect_equalities(Matrices,equalities)

% These were just found
Aeq_fixed = Matrices.G(equalities,:);
Beq_fixed = -Matrices.E(equalities,:);
beq_fixed = Matrices.W(equalities,:);

% For numerical reasons, remove fixed variables in old constraints
% (This will remove big numbers from big-M)
if ~isempty(Matrices.Aeq)
    fixed = find(Matrices.lb(1:Matrices.nu) == Matrices.ub(1:Matrices.nu));
    if ~isempty(fixed)
        Matrices.beq = Matrices.beq - Matrices.Aeq(:,fixed)*Matrices.lb(fixed);
        Matrices.Aeq(:,fixed) = 0;
    end

    skip = find(~any(full([Matrices.Aeq Matrices.beq Matrices.Beq]),2));
    if ~isempty(skip)
        Matrices.Aeq(skip,:) = [];
        Matrices.Beq(skip,:) = [];
        Matrices.beq(skip,:) = [];
    end
end

% These variables are fixed
fixed = find(Matrices.lb == Matrices.ub);
if ~isempty(fixed)
    fixed = fixed(find(fixed <= Matrices.nu));
end

Aeq = sparse(1:length(fixed),fixed,ones(length(fixed),1),length(fixed),size(Matrices.G,2));

Beq = zeros(length(fixed),size(Matrices.E,2));
beq = Matrices.lb(fixed);

% Merge everything with original equalities
Matrices.Aeq = [Aeq;Matrices.Aeq;Aeq_fixed];
Matrices.Beq = [Beq;Matrices.Beq;Beq_fixed];
Matrices.beq = [beq;Matrices.beq;beq_fixed];