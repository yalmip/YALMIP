function p = widenSmashedBox(p)
closeBad = find(abs(p.lb - p.ub) < 1e-6 & p.lb > p.ub);
if ~isempty(closeBad)
    c = (p.lb(closeBad) + p.ub(closeBad))/2;
    r = ones(length(closeBad),1)*1e-6;
    p.lb(closeBad) = c-r;
    p.ub(closeBad) = c+r;
end
