function p = addSDPextendable(p)
% A very common case is models with one signle continuous variable
% that only enter one SDP constraint. This can be exploited, as
% given an integer candidate, we can compute a feasible continuous
% variable by a quick GEVP
% To do this efficiently repeatedly, we precompute some data structures
% Waste of memory, but models are small normally
p.sdpextendable = length(p.noninteger_variables)==1 && length(p.K.s)==1 && any(p.K.s);
if p.sdpextendable
    intvars = p.integral_variables(:);
    convars = p.noninteger_variables;
    H = p.F_struc(startofSDPCone(p.K):end,:);
    H0 = reshape(H(:,1),p.K.s(1),p.K.s(1));if nnz(H0)/numel(H0)>0.5;H0 = full(H0);end
    Hx = reshape(H(:,1+convars),p.K.s(1),p.K.s(1));if nnz(Hx)/numel(Hx)>0.5;Hx = full(Hx);end
    Hz = H(:,1 + intvars);if nnz(Hz)/numel(Hz)>0.5;Hz = full(Hz);end
    % H0 + H*integers + H*continuous >= 0
    p.sdpfix.H0 = H0;
    p.sdpfix.Hx = Hx;
    p.sdpfix.Hz = Hz;
else
    p.sdpfix = [];
end