function p = addSDPextendable(p)
% A very common case is models with one single optimized continuous variable
% that only enter one SDP constraint. This can be exploited, as
% given an integer candidate, we can compute a feasible continuous
% variable by a quick GEVP
% To do this efficiently repeatedly, we precompute some data structures
% Waste of memory, but models are small normally

p.sdpextendable = 0;
if ~(length(p.K.s)==1 && any(p.K.s))    
    return
end
simplecase = length(p.noninteger_variables)==1;
if ~simplecase
    j = find(p.c);
    if nnz(p.Q)==0 && length(j)==1 && ismember(j,p.noninteger_variables)
        extendableVariable = j;
        remainingVariables = setdiff(1:length(p.c),j);
        p.sdpextendable = 1;
    end
else
    remainingVariables = p.integral_variables(:);
    extendableVariable = p.noninteger_variables;
    p.sdpextendable = 1;
end

if p.sdpextendable
    H = p.F_struc(startofSDPCone(p.K):end,:);
    H0 = reshape(H(:,1),p.K.s(1),p.K.s(1));if nnz(H0)/numel(H0)>0.5;H0 = full(H0);end
    Hx = reshape(H(:,1+extendableVariable),p.K.s(1),p.K.s(1));if nnz(Hx)/numel(Hx)>0.5;Hx = full(Hx);end
    Hz = H(:,1 + remainingVariables);if nnz(Hz)/numel(Hz)>0.5;Hz = full(Hz);end
    % H0 + H*remainingVariables + H*extendableVariable >= 0
    p.sdpfix.H0 = H0;
    p.sdpfix.Hx = Hx;
    p.sdpfix.Hz = Hz;
    p.sdpfix.convars = remainingVariables;
    p.sdpfix.forvars = extendableVariable;
else
    p.sdpfix = [];
end