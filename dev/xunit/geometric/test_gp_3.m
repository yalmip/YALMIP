function test_gp_3

sdpvar h w d

Awall  = 1;
Afloor = 1;

F = (0.5 <= h/w <= 2) + (0.5 <= d/w <= 2);
F = F + (2*(h*w+h*d) <= Awall) + (w*d <= Afloor);
F = [F, [h w] >=0];
sol = solvesdp(F,-(h*w*d))

mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal( double(h*w*d),  0.19245008957514, 1e-5);
mbg_asserttolequal( double([h w d]), [ 0.28867519677435   0.57735039354827   1.15470006291485], 1e-2);

