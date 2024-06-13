function r = fixedInfeasibleEquality(p)
% Find a == 0*x. Typically used after smashFixed
if any(p.K.f)
    r = find(~any(p.F_struc(1:p.K.f,2:end),2));
    if ~isempty(r)
        r = r(find(p.F_struc(r,1)) >= 1e-10);
    end
else
    r = [];
end