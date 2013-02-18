function test_optimizer7

X = sdpvar(2,4,2);
Y = sdpvar(2,4,2);
P = optimizer([],sum(sum(abs(X(:)-Y(:)))),[],Y,X);
Z = reshape(magic(4),2,4,2);
U = P{Z};
mbg_asserttrue( isequal(Z,U));


% Test nD parameter in cells
X = sdpvar(2,4,2);
Y1 = sdpvar(2,4,2);
Y2 = sdpvar(2,4,2);
P = optimizer([],sum(sum(abs(X(:)-Y1(:))))+sum(sum(abs(X(:)-Y2(:)))),[],{Y1,Y2},X);
Z = reshape(magic(4),2,4,2);
U = P{{Z,Z}};
mbg_asserttrue( isequal(Z,U));

% Test nD outputs in cells
X = sdpvar(2,4,2);
Y1 = sdpvar(2,4,2);
P = optimizer([],sum(sum(abs(X(:)-Y1(:)))),[],Y1,{X,2*X});
Z = reshape(magic(4),2,4,2);
U = P{Z};
mbg_asserttrue( isequal(Z,U{1}));
mbg_asserttrue( isequal(2*Z,U{2}));

