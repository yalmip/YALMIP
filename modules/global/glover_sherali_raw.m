function [L,U,infeasible] = glover_sherali_raw(a,a0,u,d)

n = length(a);
L_ = -inf(n,1);
U_ = inf(n,1);
infeasible = 0;
if u>=n || a0<=0 || all(a<=0)%any!!
    % Try reversing model
    if nargin >= 4 && d > 0
        [L_,U_,infeasible] = glover_sherali_raw(-a,-sum(a)+a0,length(a)-d);
        L = 1-U_;
        U = 1-L_;
    end
    L = L_;
    U = U_;
    return
end

% Map to sorted
[a_, loc] = sort(a,'descend');
if sum(max(a_(1:u),0)) < a0
    infeasible = 1;
    L = L_;
    U = U_;
    return
end

if any(a<0)
    u = min(u, max(find(cumsum(a_) >= a0)));
end

%From here, we only work with binary and in sorted
for j = 1:n
    s = a_(1:j);
    SN(j) = sum(s(s>0));
end
%SN = cumsum(a_);
try
    SN_u_plus_1 = sum(max(0,a_(1:u+1)));
catch
    1
end
SN_u_minus_1 = sum(max(0,a_(1:u-1)));

fixed_at_one = find(a_ > SN_u_plus_1-a0 + 1e-16);
if ~isempty(fixed_at_one)
    L_(fixed_at_one) = 1;
    U_(fixed_at_one) = 1;
end

jhat = u+1:n;
fixed_at_zero = jhat(find(a_(jhat) < a0-SN_u_minus_1 - 1e-16));
if ~isempty(fixed_at_zero)
    U_(fixed_at_zero) = 0;
    L_(fixed_at_zero) = 0;
end

% if all(a>0) && u == d
%     % Special case which can be exploited further
%     % We must set u variables true
%     SN_u_minus_1 = sum(a_(1:u-1));
%     a_ + SN_u_minus_1 
% end

L = -inf(n,1);
U = inf(n,1);
L(loc) = L_;
U(loc) = U_;