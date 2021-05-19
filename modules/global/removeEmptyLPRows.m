function p = removeEmptyLPRows(p)
if any(p.K.f) || any(p.K.l)
    % Clear complete zero rows
    zero_row = find(~any(p.F_struc,2));
    zero_row = zero_row(zero_row <= p.K.f + p.K.l);
    if ~isempty(zero_row)
        p.F_struc(zero_row,:) = [];
        p.K.l =  p.K.l - nnz(zero_row > p.K.f);
        p.K.f =  p.K.f - nnz(zero_row <= p.K.f);
    end
end