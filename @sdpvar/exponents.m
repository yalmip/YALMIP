function [ML,dummy] = exponents(poly,x)
%EXPONENTS Internal function to extract powers of nonlinear expression

mt = yalmip('monomtable');
x_lin = getvariables(poly);

x_var = [];
for i = 1:length(x)
    x_var = [x_var getvariables(extsubsref(x,i))];
end

ML = mt(x_lin,x_var);
if any(full(poly.basis(:,1)))
    ML = [zeros(1,length(x));ML];
end
dummy = [];