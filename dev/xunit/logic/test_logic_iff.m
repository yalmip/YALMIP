function test_logic_iff

yalmip('clear')

binvar x y
solvesdp(iff(x,y <= 0.5),2*x+y)
mbg_asserttolequal(double(x),0, 1e-4);
mbg_asserttolequal(double(y),1, 1e-4);
solvesdp(iff(x,y <= 0.5),x+2*y)
mbg_asserttolequal(double(x),1, 1e-4);
mbg_asserttolequal(double(y),0, 1e-4);

solvesdp(iff(x,y == 1),x+y)
mbg_asserttolequal(double(x),0, 1e-4);
mbg_asserttolequal(double(y),0, 1e-4);

solvesdp(iff(x,y == 1),-2*x+y)
mbg_asserttolequal(double(x),1, 1e-4);
mbg_asserttolequal(double(y),1, 1e-4);
solvesdp(iff(x,y == 1),x-2*y)
mbg_asserttolequal(double(x),1, 1e-4);
mbg_asserttolequal(double(y),1, 1e-4);

solvesdp(iff(x,y == 0),x)
mbg_asserttolequal(double(x),0, 1e-4);
mbg_asserttolequal(double(y),1, 1e-4);
solvesdp(iff(x,y == 0),y)
mbg_asserttolequal(double(x),1, 1e-4);
mbg_asserttolequal(double(y),0, 1e-4);

solvesdp(iff(x,y == 0),2*x+y)
mbg_asserttolequal(double(x),0, 1e-4);
mbg_asserttolequal(double(y),1, 1e-4);
solvesdp(iff(x,y == 0),x+2*y)
mbg_asserttolequal(double(x),1, 1e-4);
mbg_asserttolequal(double(y),0, 1e-4);

% Polytopic i.f.f binary equality
sdpvar x
binvar y
solvesdp([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],y,sdpsettings('debug',1))
mbg_asserttrue(abs(double(x))>.0999);
mbg_asserttolequal(abs(double(y)),0, 1e-3);
solvesdp([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],x^2,sdpsettings('debug',1))
mbg_asserttolequal(double(x),0, 1e-3);
mbg_asserttolequal(double(y),1, 1e-3);

solvesdp([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],200*x^2+y,sdpsettings('debug',1))
mbg_asserttolequal(double(x),0, 1e-3);
mbg_asserttolequal(double(y),1, 1e-3);

solvesdp([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],-x,sdpsettings('debug',1))
mbg_asserttolequal(double(x),3, 1e-3);
mbg_asserttolequal(double(y),0, 1e-3);

solvesdp([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],x,sdpsettings('debug',1))
mbg_asserttolequal(double(x),-3, 1e-3);
mbg_asserttolequal(double(y),0, 1e-3);

solvesdp([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],x-y,sdpsettings('debug',1))
mbg_asserttolequal(double(x),-3, 1e-3);
mbg_asserttolequal(double(y),0, 1e-3);


% Polytopic x equality iff binary
sdpvar x
binvar y z
solvesdp([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3<=x<=3,-10<=z<=10],y,sdpsettings('debug',1))
mbg_asserttolequal(double(y),0, 1e-3);
mbg_asserttolequal(double(z),0, 1e-3);

solvesdp([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],-y,sdpsettings('debug',1))
mbg_asserttolequal(double(y),1, 1e-3);
mbg_asserttolequal(double(z),1, 1e-3);

solvesdp([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],x,sdpsettings('debug',1))
mbg_asserttolequal(double(y),0, 1e-3);

solvesdp([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],x-y,sdpsettings('debug',1))
mbg_asserttolequal(double(y),0, 1e-3);
mbg_asserttolequal(double(x),-3, 1e-3)

solvesdp([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],x-5*y,sdpsettings('debug',1))
mbg_asserttolequal(double(y),1, 1e-3);
mbg_asserttolequal(double(z),1, 1e-3)


solvesdp([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],-y,sdpsettings('debug',1))
mbg_asserttolequal(double(y),1, 1e-3);
mbg_asserttolequal(double(z),1, 1e-3)
mbg_asserttolequal(double(x),2, 1e-3)

solvesdp([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],-z+x,sdpsettings('debug',1))
mbg_asserttolequal(double(y),0, 1e-3);
mbg_asserttolequal(double(z),1, 1e-3)
mbg_asserttolequal(double(x),-3, 1e-3)

solvesdp([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],z,sdpsettings('debug',1))
mbg_asserttolequal(double(y),0, 1e-3);
mbg_asserttolequal(double(z),0, 1e-3)

solvesdp([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],z+x,sdpsettings('debug',1))
mbg_asserttolequal(double(y),0, 1e-3);
mbg_asserttolequal(double(z),0, 1e-3)
mbg_asserttolequal(double(x),-3, 1e-3)

sdpvar z
solvesdp([iff([x <= 2, z >= 1],y == 1),-3 <= x <= 3,-10<=z<=10],z,sdpsettings('debug',1))
mbg_asserttolequal(double(y),0, 1e-3);
mbg_asserttolequal(double(z),-10, 1e-3)

solvesdp([iff([x <= 2, z >= 1],y == 1),-3 <= x <= 3,-10<=z<=10],-z+x,sdpsettings('debug',1))
mbg_asserttolequal(double(y),1, 1e-3);
mbg_asserttolequal(double(z),10, 1e-3)
mbg_asserttolequal(double(x),-3, 1e-3)



binvar d1 d2 x
solvesdp([iff(x >= 3, [d1+d2 == 2]),-5<x<5],-x,sdpsettings('debug',1))

