function q = knapsack_cover_lift_letchford(a,b,Cset)
% On lifted cover inequalities: A new lifting procedure with unusual properties
% Letchford, Souli

amagic = find_abar(a,b,Cset);
a_bar = min(a,amagic);

Nset = setdiff(1:length(a),Cset);
S = cumsum(sort(a_bar(Cset),'descend'));
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


function a_bar = find_abar(a,b,C)
l = sort(a(C),'descend');
a_bar = l(1);
c = length(C);
sigma = sum(a(C))-b;
for k = 1:c-1
    delta = a_bar  - l(k+1);
    if k*delta < sigma
        a_bar = l(k+1);
        sigma = sigma - k*delta;
    else
        a_bar = a_bar - sigma/k;
        sigma = 0;
        break
    end
end
if sigma > 0
    a_bar = b/c;
end