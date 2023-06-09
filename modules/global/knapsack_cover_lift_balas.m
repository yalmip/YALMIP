function q = knapsack_cover_lift_balas(a,Cset)
Nset = setdiff(1:length(a),Cset);
S = cumsum(sort(a(Cset),'descend'));
q = zeros(1,length(a));
q(Cset)=1;
for k = Nset
    j = max(find(a(k) >= S));
    if ~isempty(j)
        j2 = min(find(a(k) < S));
        if isempty(j2)
            'Weird, a variable should be 0'
        else
            q(k)=0;
        end
    end
end
