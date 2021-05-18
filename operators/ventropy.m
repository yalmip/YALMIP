function f = ventropy(x)

f = [];
d = size(x);
x = x(:);
for i = 1:length(x)
    f = [f;entropy(x(i))];
end
f = reshape(f,d);