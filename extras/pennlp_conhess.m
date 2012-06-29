function [nnz, row, col, val] = pennlp_fun(i,x)

try
H = datasaver(6,x,i+1);
H = reshape(H,sqrt(length(H)),sqrt(length(H)));
[row,col,val] = find(tril(H));
nnz = length(val);
row = row';
col = col';
val = val';
catch
end

%[nnz, row, col, val] = hg(i,x);

