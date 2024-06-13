function q = knapsack_cover_lift_balas(a,Cset)
Nset = setdiff(1:length(a),Cset);
S = cumsum(sort(a(Cset),'descend'));
S = [0 S(:)'];
q = zeros(1,length(a));
q(Cset)=1;
for k = Nset
    lambda = max(find(a(k) >= S))-1;
    q(k) = lambda;
end
