function [equalities,redundant] = mpt_detect_fixed_rows(Matrices)

A = [ Matrices.G -Matrices.E];
b = Matrices.W;
lb = Matrices.lb;
ub = Matrices.ub;

fixed = find(lb == ub);
if ~isempty(fixed)
    b = b-A(:,fixed)*lb(fixed);
    A(:,fixed) = [];
    ub(fixed)=[];
    lb(fixed)=[];
end

[AA, BB, AAeq, BBeq,indeq] = mpt_ineq2eq(A, b);
if ~isempty(indeq)
    i = indeq(:,1);
    j = indeq(:,2);
    equalities = unique(i);
    redundant = unique([i;j]);
    equalities = setdiff(equalities,intersect(i,j));  % The case [43 51;51 89;89 97]
    AL0A  = (A>0).*A;
    AG0A  = (A<0).*A;
    bi_up = AL0A*ub+AG0A*lb;
    redundant = unique([redundant;find([bi_up<b])]);
else
    AL0A  = (A>0).*A;
    AG0A  = (A<0).*A;
    bi_up = AL0A*ub+AG0A*lb;
    redundant = unique(find(bi_up<b));
    equalities = [];
end
