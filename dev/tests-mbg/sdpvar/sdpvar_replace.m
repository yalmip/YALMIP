function sdpvar_replace

% Checks that the 0^0 bug in MATLAB6.5 LINUX
% is avoided
sdpvar x t
p = x^2+t;
y = replace(p,t,0);
mbg_assertequal(getbase(y), getbase(x^2));
