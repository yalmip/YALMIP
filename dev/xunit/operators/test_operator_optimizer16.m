function test_operator_optimizer16

sdpvar x y
sdpvar a t
p = x^4+(x-a)^2 + a^2;
[F,objective] = compilesos(sos(p-t),-t,sdpsettings('sos.model',2),[t;a]);
P = optimizer(F,objective,sdpsettings('solver','+sdpt3'),a,t);
s1 = P{2};

a = 2;
p = x^4+(x-a)^2 + a^2;
solvesos(sos(p-t),-t)
s2 = double(t);
mbg_asserttolequal(s1,s2,1e-4);

sdpvar x
sdpvar a t
p = x^4+(x-a)^2 + a^2;
P = optimizer(sos(p-t),-t,sdpsettings('solver','+sdpt3'),a,t);
s3 = P{2};
mbg_asserttolequal(s1,s3,1e-4);
