function gamscontrol1

A = [1 2;-3 0];B = [1;1];
[K0,P0] = lqr(A,B,eye(2),1);
P = sdpvar(2,2);setsdpvar(P,2*P0);K0(K0>1)=1;K0(K0<-1)=-1;
K = sdpvar(1,2);setsdpvar(K,-K0);
F = set(K<=1)+set(K>=-1)+set(P>=0)+set((A+B*K)'*P+P*(A+B*K) <= -eye(2)-K'*K);
F = F+lmi(diag(P)>=0)+lmi(P(:)>=-151) + lmi(P(:)<=150) + lmi(P>=P0)+lmi(K>=-100) + lmi(K<=100);

obj = trace(P);

sol = solvesdp(F,obj,sdpsettings('solver','bmibnb','bmibnb.upper','penbmi,none'))

mbg_asserttolequal(sol.problem,0, 1e-5);

mbg_asserttolequal(double(obj),5.4615, 2e-2);