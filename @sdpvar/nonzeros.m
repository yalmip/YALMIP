function y = nonzeros(x)
    y = x;
    b = x.basis;
    f = find(sum(b~=0,2));
    b = b(f,:);
    y.basis = [b];
    y.dim(1) = length(f);
    y.dim(2) = 1;
end
