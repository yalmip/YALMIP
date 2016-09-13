function test_global_propagate_upper

% tTest for improved upper bound propagation, issue #319
yalmip('clear') 
x = sdpvar(2,1);
p = 100*(x(2)-x(1)^2)^2 + (1-x(1))^2;   
nlcon =  [x(1) + x(2)^2;                   
          x(1)^2 + x(2);
          x(1)^2 + x(2)^2 - 1] >= 0;
sol = optimize([-.5 <= x(1) <= .5,nlcon],p,sdpsettings('solver','bmibnb'))
mbg_asserttolequal(value(p),38.1987, 1e-2);