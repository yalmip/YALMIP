function test_sdpvar_diag

P = sdpvar(3,3,'skew');
assertTrue([0;0;0],diag(P));