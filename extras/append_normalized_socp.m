function [F_struc,K,c,Q,UB,LB,x0,Qi,Li,ri] = append_normalized_socp(F_struc,K,c,Q,UB,LB,x0)

if any(K.q)
    % We don't support initials here yet
    % x0 = [];   
    % To simplify code, we currently normalize everything to z'*z<z0^2
    % This done by writting [c^Tx+d;Ax+b] in cone as
    % [c^Tx+d;Ax+b]==[z0;z],  [z0;z] in cone
    nNew = sum(K.q);
    nOriginal = length(c);
    Ftemp = F_struc(1+K.f+K.l:end,:);
    if ~isempty(x0)
        x0 = [x0(:);Ftemp*[1;x0(:)]];
    end
    F_strucSOCP = [Ftemp -speye(nNew)];
    F_struc = [F_struc(1:K.f+K.l,:) spalloc(K.f+K.l,nNew,0)];
    F_struc = [F_strucSOCP;F_struc];
    K.f = K.f + nNew;
    c(end+nNew) = sparse(0);
    Q(end+nNew,end+nNew) = sparse(0);
    UB = [UB;inf(nNew,1)];
    lb_local = -inf(sum(K.q),1);lb_local([1 1+cumsum(K.q(1:end-1))])=0;
    LB = [LB;lb_local];
    if nargout > 7
        % Cplex interface needs explicit representations of z'Q*z+Lz+r
        top = nOriginal+1;
        ri = zeros(1,length(K.q));
        Li = spalloc(nOriginal+nNew,length(K.q),0);
        nTotal = nOriginal+nNew;
        for i = 1:length(K.q);
            ind = top:top+K.q(i)-1;
            vals = [-1 ones(1,K.q(i)-1)];
            aux = diag(sparse(ind,1,vals,nTotal,1));
            Qi{i} = aux;
            top = top + K.q(i);
        end
    end
else
    Qi = [];
    ri = [];
    Li = [];
end


