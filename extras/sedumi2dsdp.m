function [C,A,b,blk] = sedumi2dsdp(F_struc,c,K)
%SEDUMI2DSDP Internal function to convert SeDuMi structure to format needed in DSDP

nvars = size(F_struc,2)-1;
A = [];
block  = 1; 
top = 1;

if K.l>0
    blk{block,1} = 'diag';
    blk{block,2} = K.l;
    C{block,1} = sparse(F_struc(1:K.l,1));
    for var_index = 1:nvars
        A{block,var_index} = -sparse(F_struc(1:K.l,var_index+1));
    end 
    block = block+1;
    top = top+K.l;
end

% if K.q>0
%     constraints = 1;
%     while (constraints<=length(K.q))
%         n = K.q(constraints);
%         Cvec = F_struc(top:top+n-1,1);
%         C{block,1} = [Cvec(1) Cvec(2:end)';Cvec(2:end) Cvec(1)*speye(n-1,n-1)];
%         Avec = -F_struc(top:top+n-1,2:end);
%         for var_index = 1:nvars
%             lh = Avec(2:end,var_index);
%             rh = Avec(1,var_index);
%             A{block,var_index} = [rh lh';lh rh*speye(n-1)];
%         end
%         blk{block,1} = 'nondiag';
%         blk{block,2} = n;
%         constraints = constraints+1;
%         block = block+1;
%         top = top+n;
%     end
% end

if K.s>0
    constraints = 1;
    while (constraints<=length(K.s))
        n = K.s(constraints);
        Cvec = F_struc(top:top+n^2-1,1);
        C{block,1} = reshape(Cvec,n,n);
        Avec = -F_struc(top:top+n^2-1,(1:nvars)+1);
        Avec = reshape(Avec,n,n*nvars);
        left = 1;
        for var_index = 1:nvars
            A{block,var_index} = Avec(:,left:left+n-1);left = left+n;
        end
        blk{block,1} = 'nondiag';
        blk{block,2} = n;
        constraints = constraints+1;
        block = block+1;
        top = top+n*n;
    end
end
% And we solve dual...
b = -c(:);
