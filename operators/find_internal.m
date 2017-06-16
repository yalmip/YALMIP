function Model = find_internal(pos,x)

n = length(x);

% Binary to model if element is non-zero
dp = binvar(1,n); % Positve
dn = binvar(1,n); % Negative
d = binvar(1,n);

% Binary to model if element is first non-zero
e = binvar(1,n);

[M,m] = derivebounds(x);
M = max(max(abs(M),max(abs(m))));
% Activate d iff value is sufficiently non-zero
Model = [x <= 0.00001 + M*dp, x >= -0.00001 - M*dn, -M*d <= x <= M*d,
         x >= 0.00001-M*(1-dp),   x <= -0.00001 + M*(1-dn)];

% e(i) should be 1 if d(i) is the first 1
Model = [Model, d == dn + dp, e(1) == d(1), sum(e)==1];
for i = 2:n
     % No non-zeros so far so e must be 0
     Model = [Model, e(i) <= sum(d(1:i))];
     % Has to be activated if first 1
     Model = [Model, e(i) >= d(i)-sum(d(1:i-1))];
end
Model = [Model, integer(pos), pos == sum(e.*(1:n))];