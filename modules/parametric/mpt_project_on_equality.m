function [Matrices,infeasible] = mpt_project_on_equality(Matrices)
% Projects the whole mp(Q)LP problem on Aeq*U + Beq*x = beq
infeasible = 0;
if length(Matrices.beq) > 0
        
    [ii,jj,kk]=unique([Matrices.Aeq Matrices.Beq Matrices.beq],'rows');
    Matrices.Aeq = Matrices.Aeq(jj,:);
    Matrices.Beq = Matrices.Beq(jj,:);
    Matrices.beq = Matrices.beq(jj,:);
    
    if nnz(Matrices.Aeq)==0 & nnz(Matrices.Beq)>0
        % Lower dimensional in parametric space       
        infeasible = 1;
        return
    end
    
    
    
    [Qh,Rh,e] = qr(full(Matrices.Aeq),0);
    r = max(find(sum(abs(Rh),2)>1e-10));
    
   
    % The dependent
    v1 = e(1:r);
    % The basis
    v2 = e(r+1:end);

    % H1u1+H2u2 = Mv + g
    Aeq1 = Matrices.Aeq(:,v1);
    Aeq2 = Matrices.Aeq(:,v2);

    Aeqtilde = [-Aeq1\Aeq2;eye(size(Aeq2,2))];
    Beqtilde = [-Aeq1\Matrices.Beq;zeros(size(Aeq2,2),size(Matrices.Beq,2))];
    beqtilde = [Aeq1\Matrices.beq;zeros(size(Aeq2,2),1)];

    if  any(any(isnan(Aeqtilde))) | any(any(isnan(Aeqtilde))) | any(any(isnan(Aeqtilde)))
        infeasible = 1;       
        return
    end

    s = 1:size(Matrices.Aeq,2);
    p = zeros(1,length(s));
    for i = 1:length(s)
        pi = find(s(i)==e);
        if ~isempty(pi)
            p(i) = pi;
        end
    end
    % This is what we would do in ML7.1
    % [dummy,p] = ismember(1:size(Matrices.Aeq,2),e);
        
    S1 = Aeqtilde(p,:);
    S2 = Beqtilde(p,:);
    S3 = beqtilde(p,:);
    
    if  norm([Matrices.Aeq*S1 Matrices.Aeq*S2+Matrices.Beq Matrices.Aeq*S3-Matrices.beq],inf) > 1e-10
        % Lower dimensional in parametric space  
        % We basically have something like
        % u     == 1
        % u + x == 1
        infeasible = 1;
        return
    end
     
    % New parameterization U = S1*z + S2*x + S3
    M = Matrices;
    Matrices.G = M.G*S1;
    Matrices.E = M.E-M.G*S2;
    Matrices.W = M.W-M.G*S3;    
    Matrices.nu = size(Matrices.G,2);
        
    if Matrices.qp
        Matrices.H  = S1'*M.H*S1;
        Matrices.F  = M.F*S1+S2'*M.H*S1;
        Matrices.Y = M.Y + S2'*M.H*S2+0.5*(M.F*S2+S2'*M.F');
        Matrices.Cf = M.Cf*S1+S3'*M.H*S1;
        Matrices.Cc = M.Cc + M.Cf*S3;
        Matrices.Cx = M.Cx + S3'*M.F'+M.Cf*S2;
    else
        Matrices.H = M.H*S1;     
        if isfield(Matrices,'D')
            % Latest MPT does not fullify D
            Matrices.D = full((Matrices.D'*S1)');
        end
    end

    removable = find(sum(abs([Matrices.G Matrices.E Matrices.G]),2)<1e-12);
    inconsistent = intersect(removable,find(Matrices.W<-1e-10));
    if length(inconsistent)>0
        infeasible = 1;    
        return
    end

    if ~isempty(removable)
        Matrices.G(removable,:) = [];
        Matrices.E(removable,:) = [];
        Matrices.W(removable,:) = [];
    end
    
    if ~isempty(Matrices.G)
        Matrices.G(abs(Matrices.G)<eps) = 0;
    end

    % Keep the bounds for the new basis only
    Matrices.lb = [Matrices.lb(v2);Matrices.lb(end-size(Matrices.E,2)+1:end)];
    Matrices.ub = [Matrices.ub(v2);Matrices.ub(end-size(Matrices.E,2)+1:end)];

    % All equalities have been used
    Matrices.Aeq = [];
    Matrices.Beq = [];
    Matrices.beq = [];

    % This data is needed to recover original variables later
    if isempty(Matrices.getback)
        Matrices.getback.S1 = S1;
        Matrices.getback.S2 = S2;
        Matrices.getback.S3 = S3;
    else
        % This model has been reduced before, merge reductions
        oldgetback = Matrices.getback;
        Matrices.getback.S1 = oldgetback.S1*S1;
        Matrices.getback.S2 = oldgetback.S1*S2 + oldgetback.S2;
        Matrices.getback.S3 = oldgetback.S1*S3 + oldgetback.S3;
    end
end

