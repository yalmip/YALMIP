function YESNO = ismember_internal(x,p)
%ISMEMBER_INTERNAL Helper for ISMEMBER

% Author Johan Löfberg
% $Id: ismember_internal.m,v 1.1 2006-08-10 18:00:20 joloef Exp $

if isa(x,'sdpvar') & isa(p,'polytope')

    if length(p) == 1
        [H,K] = double(p);
        if min(size(x))>1
            error('first argument should be a vector');
        end
        if length(x) == size(H,2)
            x = reshape(x,length(x),1);
            YESNO = set(H*x <= K);
            return
        else
            error('Dimension mismatch.');
        end

    else
        d = binvar(length(p),1);
        YESNO = set(sum(d)==1);
        for i = 1:length(p)
            [H,K] = double(p(i));
            if min(size(x))>1
                error('first argument should be a vector');
            end
            if length(x) == size(H,2)
                x = reshape(x,length(x),1);
                lhs = H*x-K;
                [M,m] = derivebounds(lhs);
                YESNO = YESNO + set(H*x-K <= M.*(1-extsubsref(d,i)));
            else
                error('Dimension mismatch.');
            end
        end

    end
    return
end

if isa(x,'sdpvar') & isa(p,'double')

    x = reshape(x,prod(x.dim),1);
    p = p(:);

    [M,m]=derivebounds(x);

    maxp = max(p);
    minp = min(p);
    for i = 1:length(M)
        M(i) = min(M(i),maxp);
        m(i) = max(m(i),minp);
    end

    Delta = binvar(length(x),length(p),'full');
    F = set(sum(Delta,2) == 1);

    % Quick version
    xrep = repmat(x,length(p),1);
    prep = kron(p,ones(length(x),1));
    mrep = repmat(m,length(p),1);
    Mrep = repmat(M,length(p),1);
    drep = reshape(Delta,prod(size(Delta)),1);
    F = F + set((mrep-prep).*(1-drep) <= xrep-prep <= (Mrep - prep).*(1-drep));

    % Slow version
    %     for j = 1:length(x)
    %         xi = extsubsref(x,j,1);
    %        % delta = binvar(length(p),1);
    %        % F = F + set(sum(delta) == 1);
    %         for i = 1:length(p)
    %             di = extsubsref(Delta,j,i);
    %             mi = m - p(i);
    %             Mi = M - p(i);
    %             F = F + set(mi*(1-di) <= xi - p(i) <= Mi*(1-di));
    %         end
    %     end
    YESNO = F;
    return
end