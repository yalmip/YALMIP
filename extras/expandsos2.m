function  [Fout,binary_variables] = expandsos2(F,binary_variables)
sos2index = is(F,'sos2');
if ~any(sos2index)
    Fout = F;
    return
else
    LU = getbounds(F);
    Fout = F(find(~sos2index));
    for i = find(sos2index(:)')
        lambda = recover(getvariables(F(i)));
        lb = [];
        ub = [];
        for j = 1:length(lambda)
            lb(j) = LU(getvariables(lambda(j)),1);
            ub(j) = LU(getvariables(lambda(j)),2);
        end
        if any(isinf(lb))
            error('There are variables in a SOS2 constraint with no explicit lower bound');
        elseif any(isinf(ub))
            error('There are variables in a SOS2 constraint with no explicit upper bound');
        end
        n = length(lambda)-1;
        r = binvar(n,1);binary_variables = [binary_variables,getvariables(r)];
        Fout = [Fout,lb(1)*r(1) <= lambda(0+1) <= r(1)*ub(1), sum(r)==1];
        for l=1:n-1
            Fout = [Fout,lb(l+1)*(r(l)+r(l+1)) <= lambda(l+1)<= ub(l+1)*(r(l)+r(l+1))];
        end
        Fout = [Fout,lb(end)*r(end) <= lambda(end)<=r(end)*ub(end)];
    end
end