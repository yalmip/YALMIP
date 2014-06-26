function ex14_1_6

sdpvar x1 x2 x3 x4 x5 x6 x7 x8 x9 objvar; 

F = set([]);

F = F + set(  - x9 + objvar == 0);

F = F + set(0.004731*x1*x3 - 0.1238*x1 - 0.3578*x2*x3 - 0.001637*x2 - 0.9338*x4 + x7      - x9 <= 0.3571);

F = F + set( 0.1238*x1 - 0.004731*x1*x3 + 0.3578*x2*x3 + 0.001637*x2 + 0.9338*x4 - x7      - x9 <= -0.3571);

F = F + set( 0.2238*x1*x3 + 0.2638*x1 + 0.7623*x2*x3 - 0.07745*x2 - 0.6734*x4 - x7 - x9      <= 0.60220);

F = F + set( (-0.2238*x1*x3) - 0.2638*x1 - 0.7623*x2*x3 + 0.07745*x2 + 0.6734*x4 + x7      - x9 <= -0.6022);

F = F + set( x6*x8 + 0.3578*x1 + 0.004731*x2 - x9 <= 0);

F = F + set(  - x6*x8 - 0.3578*x1 - 0.004731*x2 - x9 <= 0);

F = F + set(  - 0.7623*x1 + 0.2238*x2 == -0.3461);

F = F + set( sqr(x1) + sqr(x2) - x9 <= 1);

F = F + set( (-sqr(x1)) - sqr(x2) - x9 <= -1);

F = F + set( sqr(x3) + sqr(x4) - x9 <= 1);

F = F + set( (-sqr(x3)) - sqr(x4) - x9 <= -1);

F = F + set( sqr(x5) + sqr(x6) - x9 <= 1);

F = F + set( (-sqr(x5)) - sqr(x6) - x9 <= -1);

F = F + set( sqr(x7) + sqr(x8) - x9 <= 1);

F = F + set( (-sqr(x7)) - sqr(x8) - x9 <= -1);



F = F + set( -1 <= [x1 x2 x3 x4 x5 x6 x7 x8] <= 1);

sol = solvesdp(F,objvar,sdpsettings('solver','bmibnb','bmibnb.upper','fmincon'))
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double([objvar ]),0, 1e-5);

function y = sqr(x)
y = x*x;