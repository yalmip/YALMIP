function [b,AC] = sedumi2dsdp(F_struc,c,K)
%SEDUMI2DSDP5 Internal function to convert SeDuMi structure to format needed in DSDP 5

nvars = size(F_struc,2)-1;

A = [];
block  = 1; 
top = 1;

if K.l>0    
    AC{block,1} = 'LP';
    AC{block,2} = K.l;
    A = F_struc(top:top+K.l-1,:);   
    AC{block,3} = [-A(:,2:end) A(:,1)];
    block = block+1;
    top = top+K.l;
end

if K.s>0
    for i = 1:length(K.s)
        n = K.s(i);
        AC{block,1} = 'SDP';
        AC{block,2} = n;
        A = F_struc(top:top+n^2-1,:);
        indicies = triu(reshape(1:n^2,n,n));
        indicies = indicies(find(indicies));
        AC{block,3} = [-A(indicies,2:end) A(indicies,1)];
        block = block+1;
        top = top+n*n;
    end
end
% And we solve dual...
b = -c(:);
