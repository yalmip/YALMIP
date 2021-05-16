function [p,pss] = propagate_bounds_from_qualities(p)
pss = [];
if isempty(p.nonlinear) && p.feasible
    if p.K.f>0
        Aeq = -p.F_struc(1:p.K.f,2:end);
        beq = p.F_struc(1:p.K.f,1);
        A = [Aeq;-Aeq];
        b = [beq;-beq];
        [p.lb,p.ub,redundant,pss] = tightenbounds(A,b,p.lb,p.ub,p.integer_variables,p.binary_variables,ones(length(p.lb),1));
    end
    pss=[];
    if p.K.l>0
        A = -p.F_struc(1+p.K.f:p.K.f+p.K.l,2:end);
        b = p.F_struc(1+p.K.f:p.K.f+p.K.l,1);
        [p.lb,p.ub,redundant,pss] = tightenbounds(A,b,p.lb,p.ub,p.integer_variables,p.binary_variables,ones(length(p.lb),1));
        if length(redundant)>0
            pss.AL0A(redundant,:)=[];
            pss.AG0A(redundant,:)=[];
            p.F_struc(p.K.f+redundant,:)=[];
            p.K.l = p.K.l - length(redundant);
        end
    end
end
