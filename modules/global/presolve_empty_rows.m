function p = presolve_empty_rows(p)

if p.K.f > 0
    j = find(~any(p.F_struc(1:p.K.f,:),2));
    p.F_struc(j,:)=[];
    p.K.f = p.K.f - length(j);
end

% FIXME: Add LP