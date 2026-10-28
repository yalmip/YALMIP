function bound = bnbRoundObjectiveSpectrum(spectrum,bound)
% Lift a valid lower bound to the next possible objective interval.
% A feasible cost belongs to [value-radius,value+radius]. Search with the
% lower interval endpoint of the incoming bound, then return a conservative
% lower endpoint of the first overlapping or later interval. This function
% must not be used to evaluate a feasible solution's actual objective.
if isempty(spectrum) || ~isfinite(bound)
    return
end
% Preserve the usual numerical margin on the incoming SDP bound.
margin = 1e-7*max(1,abs(bound));
target = bound-margin-spectrum.radius;
values = spectrum.values;
left = 1; right = numel(values)+1;
while left<right
    middle = floor((left+right)/2);
    if values(middle)<target
        left = middle+1;
    else
        right = middle;
    end
end
if left>numel(values)
    candidate = spectrum.cap;
else
    candidate = min(spectrum.cap,values(left)-spectrum.radius-margin);
end
bound = max(bound,candidate);
end
