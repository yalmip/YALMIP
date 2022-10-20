function [L_,U_] = glover_sherali_raw(a,a0,u)
    
L_ = -inf(length(a),1);
U_ = inf(length(a),1);
% Map to sorted
[a_, loc] = sort(a,'descend');
n = length(a_);

%From here, we only work with binary and in sorted
for j = 1:n
    s = a_(1:j);
    SN(j) = sum(s(s>0));
end
%SN = cumsum(max(a_,0));
try
    SN_u_plus_1 = sum(a_(1:u+1));
catch
    1
end
SN_u_minus_1 = sum(a_(1:u-1));

fixed_at_one = find(a_ > SN_u_plus_1-a0);
if ~isempty(fixed_at_one)
    L_(fixed_at_one) = 1;
    U_(fixed_at_one) = 1;
end

jhat = u+1:n;
fixed_at_zero = jhat(find(a_(jhat) < a0-SN_u_minus_1));
if ~isempty(fixed_at_zero)
    U_(fixed_at_zero) = 0;
    L_(fixed_at_zero) = 0;
end

L = L_(loc);
U = U_(loc);