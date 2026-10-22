function tests = test_sdpvar_review_regressions
tests = functiontests(localfunctions);

function test_prod_dimensions(testCase)
for shape = {[3 1],[1 3],[2 3]}
    dims = shape{1};
    x = sdpvar(dims(1),dims(2),'full');
    z = reshape(1:prod(dims),dims);
    assign(x,z);
    testCase.verifyEqual(value(prod(x)),prod(z));
    for dim = 1:3
        testCase.verifyEqual(value(prod(x,dim)),prod(z,dim));
    end
end

function test_prod_zero(testCase)
x = sdpvar(1);
testCase.verifyEqual(value(prod([x;0])),0);
testCase.verifyEqual(value(prod([0 x])),0);
assign(x,2);
testCase.verifyEqual(value(prod([x 0;x x])),[4 0]);

function test_median_two(testCase)
x = sdpvar(2,1);
assign(x,[2;4]);
testCase.verifyEqual(value(median(x)),3);
testCase.verifyEqual(value(median(x')),3);

function test_sort_descending(testCase)
for z = {[1 10;2 20],[1 3 2],[1;3;2]}
    a = z{1};
    x = sdpvar(size(a,1),size(a,2),'full');
    assign(x,a);
    [y,loc] = sort(x,'descend');
    [expected,indices] = sort(a,'descend');
    testCase.verifyEqual(value(y),expected);
    testCase.verifyEqual(value(loc),indices);
    testCase.verifyEqual(value(sort(x,'descend')),expected);
end

function test_trapz_matrix(testCase)
x = sdpvar(3,2,'full');
z = [1 3;2 4;7 8];
assign(x,z);
testCase.verifyEqual(value(trapz(x)),trapz(z));

function test_cumtrapz_scalar(testCase)
x = sdpvar(1);
assign(x,2);
testCase.verifyEqual(value(cumtrapz(x)),0);
testCase.verifyEqual(value(cumtrapz(3+x)),0);

function test_cumsum_scalar_and_dimensions(testCase)
x = sdpvar(1);
assign(x,3);
testCase.verifyEqual(value(cumsum(2+x)),5);
y = sdpvar(3,2,'full');
z = [1 2;3 4;5 6];
assign(y,z);
for dim = 1:4
    testCase.verifyEqual(value(cumsum(y,dim)),cumsum(z,dim));
end

function test_mrdivide_rectangular(testCase)
x = sdpvar(2,3,'full');
z = [1 2 3;4 5 6];
assign(x,z);
for divisor = {[1 0 0;0 1 0],[eye(3);1 1 1]}
    b = divisor{1};
    testCase.verifyEqual(value(x/b),z/b,'AbsTol',1e-12);
end

function test_mrdivide_infinity(testCase)
x = sdpvar(2,3,'full');
testCase.verifyEqual(x/Inf,zeros(2,3));
testCase.verifyEqual(x/(-Inf),zeros(2,3));

function test_fft_shapes_and_lengths(testCase)
x = sdpvar(3,1);
z = [1;2;4];
assign(x,z);
for row = 0:1
    if row
        a = x.'; b = z.';
    else
        a = x; b = z;
    end
    testCase.verifyEqual(value(fft(a)),fft(b),'AbsTol',1e-12);
    for n = [0 1 2 3 5]
        testCase.verifyEqual(value(fft(a,n)),fft(b,n),'AbsTol',1e-12);
    end
end

function test_jacobian_explicit_order(testCase)
x = sdpvar(3,1);
f = x(1)^2+2*x(2)^2+3*x(3)^2;
assign(x,[1;2;3]);
orders = perms(1:3);
for i = 1:size(orders,1)
    order = orders(i,:);
    expected = [2 8 18];
    testCase.verifyEqual(value(jacobian(f,x(order),1)),expected(order));
    testCase.verifyEqual(value(jacobian(f,x(order))),expected(order));
    testCase.verifyEqual(value(jacobian([f;2*f],x(order),1)),[expected(order);2*expected(order)]);
end
