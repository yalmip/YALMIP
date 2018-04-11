function test_sdpvar_subasgn

x = sdpvar(1,5);
xx = x;
x(2:3) = [1 x(1)];
t1 = isequal(x-[xx(1) 1 xx(1) xx(4:5)],zeros(1,5));

yalmip('clear')
x = sdpvar(1,5);
xx = x;
x(2:3) = [x(1) x(1)];
t2 = isequal(x-[xx(1) xx(1) xx(1) xx(4:5)],zeros(1,5));

yalmip('clear')
x = sdpvar(1,5);
xx = x;
x(2:3) = [-x(1) x(1)];
t3 = isequal(x-[xx(1) -xx(1) xx(1) xx(4:5)],zeros(1,5));


yalmip('clear')
x = sdpvar(1,5);
xx = x;
x(1:3) = [x(1)];
t4 = isequal(x-[xx(1) xx(1) xx(1) xx(4:5)],zeros(1,5));

yalmip('clear')
x = sdpvar(1,5);
xx = x;
x(1:3) = [-x(1)];
t5 = isequal(x-[-xx(1) -xx(1) -xx(1) xx(4:5)],zeros(1,5));

yalmip('clear')
x = sdpvar(1,5);
y = sdpvar(1,5);
xx = x;
x(1:5) = y;
t6 = isequal(x-y,zeros(1,5));

x = sdpvar(1,5);
y = sdpvar(1,5);
xx = x;
x(1:4) = y(1:4);
t7 = isequal(x-[y(1:4) xx(5)],zeros(1,5));


x = sdpvar(1,5);
y = sdpvar(1,5);
xx = x;
x(2:4) = y(2:4);
t8 = isequal(x-[xx(1) y(2:4) xx(5)],zeros(1,5));

x = sdpvar(2,5);
y = sdpvar(1,5);
xx = x;
x(1,2:4) = y(2:4);
t9 = isequal(x-[xx(1,1) y(2:4) xx(1,5);xx(2,:)],zeros(2,5));

x = sdpvar(2,5);
y = sdpvar(1,5);
xx = x;
x(1,2:4) = 0;
t10 = isequal(x-[xx(1,1) 0 0 0 xx(1,5);xx(2,:)],zeros(2,5));

x = sdpvar(2,5);
y = sdpvar(1,5);
xx = x;
x(1:2,2) = 0;
t11 = isequal(x-[xx(:,1) [0;0] xx(:,3:end)],zeros(2,5));

x = sdpvar(2,5);
y = sdpvar(1,5);
xx = x;
x(1:2,2:end) = 0;
t12 = isequal(x-[xx(:,1) zeros(2,4)],zeros(2,5));


x = sdpvar(2,5);
y = sdpvar(1);
xx = x;
x(1,7) = y;
t13 = isequal(x-[xx [0 y;0 0]],zeros(2,7));

% Bug 1006
p = sdpvar(3,1);
p([1 1])=[4 5];
t14 = isa(p(1),'double') & isequal(p(1),5);

% Bug #190
x = sdpvar(2);
y = x;
y(:,2) = [];
t15 = isequal(y-x(:,1),[0;0]);

assertTrue(t1 & t2 & t3 & t4 & t5 & t6 & t7 & t8 & t9 & t10 &t11 & t12 & t13 & t14 & t15);