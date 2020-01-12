function tests = test_operator_optimizer6
tests = functiontests(localfunctions);

function test1(dummy)

X = sdpvar(2,3,'full','complex');
Y = sdpvar(2,3,'full','complex');
P = optimizer([],sum(sum(abs(real(X)-real(Y)))) + sum(sum(abs(imag(X)-imag(Y)))),[],Y,X);

U = reshape(1:6,2,3);
Z = P{U};
assert(norm(U-Z)<= 1e-10);

X = sdpvar(2,3,2,'full','complex');
Y = sdpvar(2,3,'full','complex');
obj = sum(sum(abs(real(X(:,:,1))-real(Y)))) + sum(sum(abs(imag(X(:,:,1))-imag(Y))));
obj = obj + sum(sum(abs(imag(X(:,:,2))-2*real(Y)))) + sum(sum(abs(real(X(:,:,2))-2*imag(Y))));
P = optimizer([],obj,[],Y,X);
U = reshape(1:6,2,3);
Z = P{U};
assert(norm(U-Z(:,:,1))<1e-10);
assert(norm(2*(U)-imag(Z(:,:,2)))<1e-10);

U = reshape(1:6,2,3) + sqrt(-1)*reshape(7:12,2,3);
Z = P{U};
assert(norm(U-Z(:,:,1)) <= 1e-5);
assert(norm(2*real(U)-imag(Z(:,:,2))) <= 1e-5);
assert(norm(2*imag(U)-real(Z(:,:,2))) <= 1e-5);

X = sdpvar(2,3,'full','complex');
Y = sdpvar(2,3,'full');
P = optimizer([],sum(sum(abs(real(X)-Y))) + 2*sum(sum(abs(imag(X)-Y))),[],X,Y);
U = reshape(1:6,2,3) + sqrt(-1)*reshape(7:12,2,3);
Z = P{U};
assert(norm(Z(:,:,1)-reshape(7:12,2,3))<1e-7);

X1 = sdpvar(2,3,'full','complex');
X2 = sdpvar(2,3,'full');
Y = sdpvar(2,3,'full');
P = optimizer([],sum(sum(abs(real(X1)-Y)))+2*sum(sum(abs(imag(X1)-Y)))+sum(sum(X2-Y)),[],{X1,X2},Y);
U1 = reshape(1:6,2,3) + sqrt(-1)*reshape(7:12,2,3);
U2 = -reshape(1:6,2,3);
Z = P{{U1,U2}};
assert(norm(Z-reshape(7:12,2,3))<1e-8);
