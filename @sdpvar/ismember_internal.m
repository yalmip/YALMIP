function YESNO = ismember_internal(x,p)
%ISMEMBER_INTERNAL Helper for ISMEMBER

% Author Johan Löfberg
% $Id: ismember_internal.m,v 1.10 2007-08-08 08:14:28 joloef Exp $

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
            disp('The polytope in the ismember condition has wrong dimension')
            error('Dimension mismatch.');
        end

    else
        d = binvar(length(p),1);
        YESNO = set(sum(d)==1);
        [temp,L,U] = bounding_box(p(1));
        for i = 1:length(p)
            [temp,Li,Ui] = bounding_box(p(i));
            L = min([L Li],[],2);
            U = max([U Ui],[],2);
        end
        for i = 1:length(p)
            [H,K] = double(p(i));
            if min(size(x))>1
                error('first argument should be a vector');
            end
            if length(x) == size(H,2)
                x = reshape(x,length(x),1);
                lhs = H*x-K;
                % Derive bounds based on YALMIPs knowledge on bounds on
                % involved variables
                [M,m] = derivebounds(lhs);
                % Strengthen by using MPTs bounding box
                %[temp,L,U] = bounding_box(p(i));
                Hpos = (H>0);
                Hneg = (H<0);
                M = min([M (H.*Hpos*U+H.*Hneg*L-K)],[],2);
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
   
    
    if numel(p)==1
        F = set(x == p);
    else
        if size(p,1)==length(x) & size(p,2)>1
            Delta = binvar(size(p,2),1);
            F = [sum(Delta) == 1, x == p*Delta];
        else
            p = p(:);
            Delta = binvar(length(x),length(p),'full');
            F = [sum(Delta,2) == 1, x == Delta*p];
        end
    end

    YESNO = F;
    return
end