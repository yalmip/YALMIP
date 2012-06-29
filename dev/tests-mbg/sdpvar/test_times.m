function test_times

% bug reported by Erik Johanesson.  Missing complex transpose in times
x = sdpvar(2,1,'full','complex');
w = randn(2,1) + sqrt(-1)*randn(2,1);
v = randn(2,1) + sqrt(-1)*randn(2,1);

assign(x,w);

double(w.*v - x.*v)

mbg_asserttrue(all(double(w.*v - x.*v) < 1e-12));