function [xtemp,fail] = sdpextendsolution(p,xtemp)
% Given z and H0 + H*x + Hz*z try to find scalar x
% such that SDP is feasible
% Used in MISDP heuristics
fail = 1;
if p.sdpextendable  
    Hy = p.sdpfix.H0 + reshape(p.sdpfix.Hz*xtemp(p.integral_variables),p.K.s(1),p.K.s(1));
    s = eig(full(p.sdpfix.Hx),full(Hy));
    s(isinf(s))=[];
    s(isnan(s))=[];
    if any(s)
        xtemp(p.noninteger_variables) = min(-1./s(s~=0));
        fail = 0;
    end
end