function test_kronecker

a = [1 2 3;4 5 6];
b = [7 8 9];

A = sdpvar(2,3);
B = sdpvar(1,3);
assign(A,a);
assign(B,b);

% We have special code for vector arguments.
% Make sure these work...
% Row vector
assertEqual(double(kron(A,b)),kron(a,b))
assertEqual(double(kron(b,A)),kron(b,a))
assertEqual(double(kron(a,B)),kron(a,b))
assertEqual(double(kron(B,a)),kron(b,a))
assertEqual(double(kron(A,B)),kron(a,b))
% Column vector
b = [7 8 9]';
B = sdpvar(3,1);
assign(B,b);
assertEqual(double(kron(A,b)),kron(a,b))
assertEqual(double(kron(b,A)),kron(b,a))
assertEqual(double(kron(a,B)),kron(a,b))
assertEqual(double(kron(B,a)),kron(b,a))
assertEqual(double(kron(A,B)),kron(a,b))

% Scalar
b = 1;
B = sdpvar(1,1);
assign(B,b);
assertEqual(double(kron(A,b)),kron(a,b))
assertEqual(double(kron(b,A)),kron(b,a))
assertEqual(double(kron(a,B)),kron(a,b))
assertEqual(double(kron(B,a)),kron(b,a))
assertEqual(double(kron(A,B)),kron(a,b))

% General matrix
b = [6 7 8;9 10 11];
B = sdpvar(2,3);
assign(B,b);
assertEqual(double(kron(A,b)),kron(a,b))
assertEqual(double(kron(b,A)),kron(b,a))
assertEqual(double(kron(a,B)),kron(a,b))
assertEqual(double(kron(B,a)),kron(b,a))
assertEqual(double(kron(A,B)),kron(a,b))


