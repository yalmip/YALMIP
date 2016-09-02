function sys=relaxvalue(X)
%RELAXVALUE Return numerial value treating nonlinear variables as independent

solution = yalmip('getsolution');
lmi_variables = X.lmi_variables;
opt_variables = solution.variables;

values = zeros(1+length(lmi_variables),1);
values(1)=1;
for i=1:length(lmi_variables)
    opt_index = find(lmi_variables(i)==opt_variables);
    if isempty(opt_index)
        values(i+1,1)=NaN;
    else
        values(i+1,1)=solution.optvar(opt_index);
    end
end

sys = X.basis*values;
sys = full(reshape(sys,X.dim(1),X.dim(2)));