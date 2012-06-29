function sdpvar_diag

% Bugs reported from Stefano di Cairano
P = sdpvar(3,3,'skew');
mbg_assertequal([0;0;0],diag(P));