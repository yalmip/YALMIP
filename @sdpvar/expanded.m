function X = expanded(Y,state)

if nargin == 1
    X = Y.extra.expanded;
else
    X = Y;
    if state == 0
        X.extra.expanded = [];
    else
        X.extra.expanded = X.lmi_variables;
    end
end
