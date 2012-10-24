function test_hessian

sdpvar x y
mbg_asserttrue(all(all(zeros(2)==hessian(x,[x y]))));
mbg_asserttrue(all(all(zeros(2)==hessian(x+y,[x y]))));
mbg_asserttrue(all(all([0 1;1 0]==hessian(x*y,[x y]))));
mbg_asserttrue(all(all([2 0;0 0]==hessian(x^2,[x y]))));
mbg_asserttrue(all(all([2]==hessian(x^2+y^2,[x]))));
mbg_asserttrue(all(all([2]==hessian(x^2+y^2,[y]))));
mbg_asserttrue(all(all([0]==hessian(y^2,[x]))));






 
