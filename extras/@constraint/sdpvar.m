function F = sdpvar(X)
% Internal class for constraint list

F = sdpvar(lmi(X));
