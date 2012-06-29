function dfdx_sort = map_to_original(dfdx,x, x_indep);

for i = 1:length(x(:))
    x_origvar(i) = getvariables(x(i));
end
x_var = getvariables(x_indep);
dfdx_sort = [];
for i = 1:length(x_origvar)
    dfdx_sort = [dfdx_sort dfdx(:,find(x_origvar(i) == x_var))];
end