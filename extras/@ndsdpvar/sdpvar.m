function F = sdpvar(X)
% SDPVAR Converts a NDSDPVAR variable to standard SDPVAR

F = sdpvar(prod(X.dim),1,[],X.lmi_variables,X.basis);