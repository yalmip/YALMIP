function x = round_solution(x,p)

x(p.integer_variables) = round(x(p.integer_variables));
x(p.binary_variables) = round(x(p.binary_variables));

