function cut = knapsack_cheapcut(row)
% Dedicated easy cut for sum ai xi >= 0
% which appears effective on some models
if all(row(2:end)>=0) & row(1)<=0
    a = row(2:end);
    b = -row(1);
    [val,loc] = sort(a,'ascend');
    g = cumsum(val);
    m = max(find(g < -row(1)));
    ok = setdiff(loc,loc(1:m));
    b_cut = -1;
    a_cut = spalloc(1,length(val),0);
    a_cut(ok) = 1;
    cut = [b_cut a_cut];
else
    cut = [];
end