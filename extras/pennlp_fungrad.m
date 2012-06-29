function  [nnz,ind,val] = pennlp_fungrad(x)


f = datasaver(2,x);
[ind,dummy,val] = find(f(:));
nnz = length(ind);