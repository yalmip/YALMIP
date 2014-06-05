function [mDIM,nBLOCK,bLOCKsTRUCT,c,F] = sedumi2sdpa(F_struc,c,K);
%SEDUMI2SDPA Internal function to convert SeDuMi structure to format needed in SDPA

start = 1;

% This is a hack. K.f is only available when called from bnb with dynamically added equalities
if K.f>0
    F_struc = [-F_struc(1:K.f,:);F_struc];
    K.l = K.l + K.f;
    K.f = 0;
end
if K.s>0
    nc = length(K.s);
else
    nc=0;
end
if K.l>0
    nc = nc+1;
end

F = cell(nc,size(F_struc,2));

% Linear constraints
bLOCKsTRUCT =[];
nl = 0;
if K.l>0
    F{1,1}=sparse(-F_struc(1:K.l,1));
    for j = 2:size(F_struc,2)      
        F{1,j}=sparse(F_struc(1:K.l,j)); 
    end
    start = K.l+1;
    bLOCKsTRUCT = [-K.l];
    nl = 1;
end

% Semidefinite constraints
if K.s>0
    for i = 1:length(K.s)
        z = tril(ones(K.s(i)),-1);z = find(z(:));
        theend = start+power(K.s(i),2)-1;
        temp = sparse(F_struc(start:theend,:))';
        temp(:,z) = 0;
        temp = temp';        
        %F{i+nl,1}=triu(sparse(-reshape(F_struc(start:theend,1),K.s(i),K.s(i))));
        F{i+nl,1}=((-reshape(temp(:,1),K.s(i),K.s(i))));
        FastZero = spalloc(K.s(i),K.s(i),0);
        nonZeros = setdiff(find(any(temp,1)),1);
        for j = nonZeros%2:size(F_struc,2)
            %bbb=reshape(F_struc(start:theend,j),K.s(i),K.s(i));
            r = temp(:,j);           
            bbb=reshape(r,K.s(i),K.s(i));
            aaa = (bbb);
            F{i+nl,j}=(aaa);           
        end
       % Zeros = setdiff(find(~any(temp,1)),1);
       
       % for j = Zeros
       %     F{i+nl,j} = [];%FastZero;            
       % end
        
        start = theend+1;
        bLOCKsTRUCT = [bLOCKsTRUCT K.s(i)];
    end
end
mDIM = size(F_struc,2)-1;
nBLOCK = length(bLOCKsTRUCT);