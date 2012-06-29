function sys=relaxdouble(X)
%RELAXDOUBLE Return numerial value treating nonlinear variables as independent

% Author Johan Löfberg 
% $Id: relaxdouble.m,v 1.1 2006-08-10 18:00:22 joloef Exp $  

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