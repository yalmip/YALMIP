function p = clean_bounds(p)
close = 1e-10>abs(p.ub - round(p.ub));
p.ub(close) = round(p.ub(close));

close = 1e-10>abs(p.lb - round(p.lb));
p.lb(close) = round(p.lb(close));
p.ub(p.binary_variables) = floor(p.ub(p.binary_variables) + 1e-2);
p.lb(p.lb<-1e12) = -inf;
p.ub(p.ub>1e12) = inf;