function Model = isoutside_internal(P)
%ISOUSIDE_INTERNAL Helper for ISOUTSIDE

% Get representation A*x <= b
B = getbase(sdpvar(P));
b = B(:,1);
A = -B(:,2:end);
x = recover(P);

% Derive big-M for some Ai*x-bi >= 0
d = binvar(length(b),1);
violation = A*x-b;
[M,m] = derivebounds(violation);
Model = [sum(d)==1, violation >= m.*(1-d)];
