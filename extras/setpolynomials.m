function setpolynomials
%SETPOLYNOMIALS Internal function

solution = sdpvar('getSolution');
lmi_variables = solution.variables;
opt_variables = solution.optvar;
sqrList = yalmip('sqrvariables');
polyvals = [];

if ~(isempty(lmi_variables) | isempty(sqrList))
    for i = 1:length(sqrList)
        indx1 = find(lmi_variables == sqrList(i,2));
        indx2 = find(lmi_variables == sqrList(i,3));
        if  ~isempty(indx1) & ~isempty(indx2)
            polyvals = [polyvals; sqrList(i,1) opt_variables(indx1)*opt_variables(indx2)];
        end    
    end
    if ~isempty(polyvals)
        assign(recover(polyvals(:,1)),polyvals(:,2));
    end
end
