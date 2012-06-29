function test_logic_mixed

binvar d1 d2

sol =  solvesdp(~(d1==d2),d2)
mbg_asserttrue(sol.problem == 0);
mbg_asserttolequal(double(d1 + d2),1,1e-3)

sol =  solvesdp((d1~=d2),d1+d2)
mbg_asserttrue(sol.problem == 0);
mbg_asserttolequal(double(d1 + d2),1,1e-3)

sol =  solvesdp(iff(d1,d2),d1+d2)
mbg_asserttrue(sol.problem == 0);
mbg_asserttolequal(double(d1 + d2),0,1e-3)

sol =  solvesdp(iff(d1,d2),d1)
mbg_asserttrue(sol.problem == 0);
mbg_asserttolequal(double(d2),0,1e-3)