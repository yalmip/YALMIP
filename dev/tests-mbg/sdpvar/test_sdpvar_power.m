function sdpvar_replace

% Bugs reported from Didier Henrion
x = sdpvar(1);
y = ([1 1 1]*x).^(0:2);
mbg_assertequal(getbase(y), getbase([1 x x^2]));

x = sdpvar(1);
y = x.^(0:2);
mbg_assertequal(getbase(y), getbase([1 x x^2]));

