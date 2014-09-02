function test_logic_sat_1

binvar aa ba ca da ea fa % NOTE a b c ... generates error when run as function !
F = (true((aa & ba & ca) | (da & ea & fa)));
sol = solvesdp(F);
mbg_asserttolequal(sol.problem,0);
mbg_asserttrue( all(double([aa ba ca])== [1 1 1])  | all(double([da ea fa])== [1 1 1]) )

