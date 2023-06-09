function [F,M,m] = max_integer_model(X,t,M,m)

if nargin < 3
    [M,m] = derivebounds(X);
end

if all(M==m) 
    F = [t == max(M)];
    return
end

if length(X) > 10
    % Integer model of max(X) has O(n^2) complexity in memory use as there
    % are n constrants with n non-zeros. Recursive variant has O(nlog(n))
    mid = floor(length(X)/2);
    t1 = sdpvar(1);
    t2 = sdpvar(1);
    % t1 = max of first half
    F1 = max_integer_model(X(1:mid),t1,M(1:mid),m(1:mid));
    % t2 = max of second half
    F2 = max_integer_model(X(mid+1:end),t2,M(mid+1:end),m(mid+1:end));
    % and take max of those to max, t = max(t1,t2)
    M3 = [max(M(1:mid));max(M(mid+1:end))];
    m3 = [min(m(1:mid));min(m(mid+1:end))];
    F3 = max_integer_model([t1;t2],t,M3,m3);
    F = F1 + F2 + F3;   
    return
end

n = length(X);
d = binvar(n,1);
F = (sum(d)==1);
F = F + (-(max(M)-min(m))*(1-d) <= t-X <= (max(M)-min(m))*(1-d));
kk = [];
ii = [];
for i = 1:n
    k = [1:1:i-1 i+1:1:n]';
    ii = [ii;repmat(i,n-1,1)];
    kk = [kk;k];
    Mm = M(k)-m(i);
end
xii = extsubsref(X,ii);
dii = extsubsref(d,ii);
xkk = extsubsref(X,kk);
F = F + (xkk <= xii+(M(kk)-m(ii)).*(1-dii));