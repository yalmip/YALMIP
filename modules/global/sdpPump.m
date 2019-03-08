function [x,upper,possible,p,successful] = sdpPump(p,x,tolerance)

upper = inf;
intvars = [p.integer_variables(:);p.binary_variables(:)];
convars = p.noninteger_variables;
possible = 0;
successful = 0;

if length(convars)==1 && length(p.K.s)==1 && p.K.s(1)>0 && sum(p.K.q)==0
    possible = 1;
    if isempty(p.sdpPumpData)
        p.sdpPumpData.H = p.F_struc(1+p.K.f+p.K.l:end,:);
        p.sdpPumpData.H0 = reshape(p.sdpPumpData.H(:,1),p.K.s(1),p.K.s(1));
        p.sdpPumpData.Hx = reshape(p.sdpPumpData.H(:,1+convars),p.K.s(1),p.K.s(1));
        p.sdpPumpData.Hz = p.sdpPumpData.H(:,1 + intvars);if nnz(p.sdpPumpData.Hz)/numel(p.sdpPumpData.Hz)>0.5;p.sdpPumpData.Hz = full(p.sdpPumpData.Hz);end    
    end
    Hy = p.sdpPumpData.H0 + reshape(p.sdpPumpData.Hz*x(intvars),p.K.s(1),p.K.s(1));
    s = eig(full(p.sdpPumpData.Hx),full(Hy));
    s(isinf(s))=[];
    s(isnan(s))=[];
    if any(s)
        x(convars) = min(-1./s(s~=0));
        if ~isnan(x(convars))
            x(convars) = max(x(convars),p.lb(convars));
            x(convars) = min(x(convars),p.ub(convars));
            if checkfeasiblefast(p,x,tolerance)
                upper = computecost(p.f,p.c,p.Q,x,p);
                successful = 1;
            else
                upper = inf;
                successful = 0;
            end
        end
    end
end
    
