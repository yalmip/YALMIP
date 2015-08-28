function info = lmiinfo(F)

info.sdp = [];
info.lin = [];
info.equ = [];
info.soc = [];
info.rlc = [];
info.pow = [];
F = flatten(F);
Counter = length(F.LMIid);
sdptop = 1;
info.sdp = zeros(Counter,3);
for i = 1:Counter
    Fi = F.clauses{i}.data;
    switch  F.clauses{i}.type;
        case {1,9,40}
            % Slightly faster, quick fix for user case...
            info.sdp(sdptop,:) = [size(Fi,1) size(Fi,2) F.LMIid(i)];
            sdptop = sdptop  + 1;
        case 2
            info.lin = [info.lin;size(Fi,1) size(Fi,2) F.LMIid(i)];
        case 3
            info.equ = [info.equ;size(Fi,1) size(Fi,2) F.LMIid(i)];
        case {4,54}
            info.soc = [info.soc;size(Fi,1) size(Fi,2) F.LMIid(i)];
        case 5
            info.rlc = [info.rlc;size(Fi,1) size(Fi,2) F.LMIid(i)];
        case 20
            info.pow = [info.pow;size(Fi,1) size(Fi,2) F.LMIid(i)];
        case {7,8}
        otherwise
            error('Error in lmiinfo. Please report bug');
    end
end
info.sdp = info.sdp(1:sdptop-1,:);