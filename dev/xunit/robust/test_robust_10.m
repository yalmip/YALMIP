function test_robust_10
% Test of Bertsime-Sim D-norm
yalmip('clear')
p = 1.5;
a = randn(2,1);

% Analytic value
o1 = max(norm(a,inf),norm(a,1)/p);

% Robust optimization over explicitly reduced model
w = sdpvar(2,1);
sdpvar t
W1 = [norm(w,1)<=1];
W2 = [norm(w,inf)<=1/p];
H = hull(W1,W2);
W = projection(H,w);
solvesdp([uncertain(w),W,a'*w<=t],t);
o2 = double(t);

% Robust optimization over original model
% FAILS, issues #123
if 0
yalmip('clear')
w = sdpvar(2,1);
sdpvar t
W1 = [norm(w,1)<=1];
W2 = [norm(w,inf)<=1/p];
[H,tt,y] = hull(W1,W2);
solvesdp([uncertain(w),H,a'*w<=t],t)
o3 = double(t);
mbg_asserttolequal(o1-o3,0, 1e-4);


% Robust optimization over original model
yalmip('clear')
w = sdpvar(2,1);
sdpvar t
W1 = [norm(w,1)<=1];
W2 = [norm(w,inf)<=1/p];
solvesdp([uncertain(w),[W1,W2],a'*w<=t],t)

o3 = double(t);
mbg_asserttolequal(o1-o3,0, 1e-4);
end


