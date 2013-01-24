function [F_struc,K,c,Q,UB,LB,Qi,Li,ri] = append_normalized_socp(F_struc,K,c,Q,UB,LB)

if K.q(1)>0
    % To simplify code, we currently normalize everything to z'*z<z0^2
    % This done by writting [c^Tx+d;Ax+b] in cone as
    % [c^Tx+d;Ax+b]==[z0;z],  [z0;z] in cone
    nNEW = sum(K.q);
    nORIGINAL = length(c);
    Ftemp = F_struc(1+K.f+K.l:end,:);
    F_strucSOCP = [Ftemp -speye(nNEW)];
    F_struc = [F_struc(1:K.f+K.l,:) spalloc(K.f+K.l,nNEW,0)];
    F_struc = [F_strucSOCP;F_struc];
    K.f = K.f + nNEW;
    c = [c;spalloc(nNEW,1,0)];
    Q = blkdiag(Q,spalloc(nNEW,nNEW,0));
    UB = [UB;inf(nNEW,1)];
    for i = 1:length(K.q);
        LB = [LB;0;-inf(K.q(i)-1,1)];
    end
    if nargout > 6
        % Cplex interface needs explicit representations of z'Q*z+Lz+r
        iCone = nORIGINAL+1;
        ri = zeros(1,length(K.q));
        Li = spalloc(nORIGINAL+nNEW,length(K.q),0);
        nTOT = nORIGINAL+nNEW;
        for i = 1:length(K.q);
            ind = iCone:iCone+K.q(i)-1;
            vals = [-1 ones(1,K.q(i)-1)];
            aux = diag(sparse(ind,1,vals,nTOT,1));
            Qi{i} = aux;
            iCone = iCone + K.q(i);
        end
    end
else
    Qi = [];
    ri = [];
    Li = [];
end


