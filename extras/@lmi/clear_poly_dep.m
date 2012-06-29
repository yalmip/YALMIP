function F = clear_poly_dep(F,x,order)

F.clauses{1}.data = clear_poly_dep(F.clauses{1}.data,x,order);