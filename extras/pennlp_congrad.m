function [nnz,ind, val] = pennlp_congrad(i,x)


try
    f = datasaver(5,x,i+1);
    [ind,dummy,val] = find(f(:));
    nnz = length(ind);
catch
  
end

%[nnz,ind, val] = dg(i,x); 

