function upper_bound = upper_from_sdpcuts(p,sdpCuts)
% In models here we have a simple unit continuous objective
% we can use sdp cuts collected on the way to quickly derive
% a crude bound

% Fix given variables
S = sdpCuts.F_struc;
k = find(p.ub == p.lb);
s = S(:,1+k)*p.ub(k);
S(:,1) = S(:,1) + s;
S(:,1+k)=[];

% Use cardinality constraints to strengthen bound
n_fixed = nnz(p.lb==1);
n_free = p.globalcardinality.up - n_fixed;
upper_bound = -inf;
for i = 1:size(S,1)
    if n_free == 0
        upper_bound = max(upper_bound,(S(i,1))./S(i,end));
    else
        row = max(S(i,2:end-1),0);
        upper_bound = max(upper_bound,(S(i,1) + sumk(row,n_free))./S(i,end));
    end
end