function setpolynomials
%SETPOLYNOMIALS Internal function

% Author Johan Löfberg
% $Id: setpolynomials.m,v 1.3 2005-04-29 08:05:01 joloef Exp $

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
