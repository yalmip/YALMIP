function test_optimizer6

X = sdpvar(2,3,'full','complex');
Y = sdpvar(2,3,'full','complex');
P = optimizer([],sum(sum(abs(real(X)-real(Y)))) + sum(sum(abs(imag(X)-imag(Y)))),[],Y,X);

U = reshape(1:6,2,3);
Z = P{U};
mbg_asserttrue(isequal(U,Z));

U = reshape(1:6,2,3) + sqrt(-1)*reshape(7:12,2,3);
Z = P{U};
mbg_asserttrue(isequal(U,Z));

X = sdpvar(2,3,'full','complex');
Y = sdpvar(2,3,'full');
P = optimizer([],sum(sum(abs(real(X)-Y))) + 2*sum(sum(abs(imag(X)-Y))),[],X,Y);
U = reshape(1:6,2,3) + sqrt(-1)*reshape(7:12,2,3);
Z = P{U};
mbg_asserttrue(isequal(Z,reshape(7:12,2,3)));

X1 = sdpvar(2,3,'full','complex');
X2 = sdpvar(2,3,'full');
Y = sdpvar(2,3,'full');
P = optimizer([],sum(sum(abs(real(X1)-Y)))+2*sum(sum(abs(imag(X1)-Y)))+sum(sum(X2-Y)),[],{X1,X2},Y);
U1 = reshape(1:6,2,3) + sqrt(-1)*reshape(7:12,2,3);
U2 = -reshape(1:6,2,3);
Z = P{{U1,U2}};
mbg_asserttrue(isequal(Z,reshape(7:12,2,3)));
