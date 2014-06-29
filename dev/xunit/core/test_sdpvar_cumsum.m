function test_sdpvar_cumsum

dx = sdpvar(1,2);
cs1 = [0 cumsum(dx) 1];
cs2 = [0 dx(1) dx(1)+dx(2) 1];
assertTrue(isequal(cs1-cs2,[0 0 0 0]))
