function tests = test_logic_sat1
tests = functiontests(localfunctions);

function test1(dummy)

binvar aa ba ca da ea fa % NOTE a b c ... generates error when run as function !
F = (true((aa & ba & ca) | (da & ea & fa)));
sol = solvesdp(F);
assert(sol.problem == 0);
assert( all(double([aa ba ca])== [1 1 1])  | all(double([da ea fa])== [1 1 1]) )

