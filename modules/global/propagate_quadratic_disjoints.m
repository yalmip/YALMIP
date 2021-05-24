function p = propagate_quadratic_disjoints(p,spliton)

if any(p.quadraticdisjoints)
   loc = find(spliton == p.quadraticdisjoints(1,:));
   if any(loc)
       % x >= u or x <= l
       u = p.quadraticdisjoints(2,loc);
       l = p.quadraticdisjoints(3,loc);
       if p.lb(spliton) > l
           p.lb(spliton) = max(p.lb(spliton),u);
       end
       if p.ub(spliton) < u
           p.ub(spliton) = min(p.ub(spliton),l);
       end
   end
end