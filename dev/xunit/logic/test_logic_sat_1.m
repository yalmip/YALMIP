function tests = test_logic_sat1
tests = functiontests(localfunctions);

function test1(dummy)

binvar aa ba ca da ea fa % NOTE a b c ... generates error when run as function !
F = (true((aa & ba & ca) | (da & ea & fa)));
sol = optimize(F);
assert(sol.problem == 0);
assert( all(value([aa ba ca])== [1 1 1])  | all(value([da ea fa])== [1 1 1]) )

