function [F_struc,K,c,Q,UB,LB] = append_normalized_socp(F_struc,K,c,Q,UB,LB)

if K.q(1)>0
    n_original = length(c);
    
    
    % To simplify code, we currently normalize everything to z'*z<x0^2
    nNEW = sum(K.q);
    
    Ftemp = F_struc(1+K.f+K.l:end,:);
    F_strucSOCP = [Ftemp -speye(nNEW)];
    F_struc = [F_struc(1:K.f+K.l,:) spalloc(K.f+K.l,nNEW,0)];
    UB = [UB;inf(nNEW,1)];
    c = [c;spalloc(nNEW,1,0)];
    Q = blkdiag(Q,spalloc(nNEW,nNEW,0));
    
    iCone = n_original+1;
    ri = zeros(1,length(K.q));
    Li = spalloc(n_original+nNEW,length(K.q),0);
    for i = 1:length(K.q);
        Qi{i} = sparse(iCone:iCone+K.q(i)-1,iCone:iCone+K.q(i)-1,[-1 ones(1,K.q(i)-1)],n_original+nNEW,n_original+nNEW);
        LB = [LB;0;-inf(K.q(i)-1,1)];
        iCone = iCone + K.q(i);
    end
    F_struc = [F_strucSOCP;F_struc];
    K.f = K.f + nNEW;
    
else
    Qi = [];
    ri = [];
    Li = [];
end


