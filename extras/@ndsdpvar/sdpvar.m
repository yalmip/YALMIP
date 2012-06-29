function F = sdpvar(X)
% SDPVAR Converts a NDSDPVAR variable to standard SDPVAR

% Author Johan Löfberg
% $Id: sdpvar.m,v 1.3 2006-07-13 19:40:59 joloef Exp $

F = sdpvar(prod(X.dim),1,[],X.lmi_variables,X.basis);