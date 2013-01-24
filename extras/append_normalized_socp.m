function [F_struc,K,c,Q,UB,LB] = append_normalized_socp(F_struc,K,c,Q,UB,LB)

if K.q(1)>0
    % To simplify code, we currently normalize everything to z'*z<x0^2
    nNEW = sum(K.q);
    
    Ftemp = F_struc(1+K.f+K.l:end,:);
    F_strucSOCP = [Ftemp -speye(nNEW)];
    F_struc = [F_struc(1:K.f+K.l,:) spalloc(K.f+K.l,nNEW,0)];
    UB = [UB;inf(nNEW,1)];
    c = [c;spalloc(nNEW,1,0)];
    Q = blkdiag(Q,spalloc(nNEW,nNEW,0));
    
    ri = zeros(1,length(K.q));
    for i = 1:length(K.q);   
        LB = [LB;0;-inf(K.q(i)-1,1)];
    end
    F_struc = [F_strucSOCP;F_struc];
    K.f = K.f + nNEW;    
end


