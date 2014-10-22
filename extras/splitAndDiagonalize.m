function NewConstraints = splitAndDiagonalize(constraints,nmin,ratio)

kept = ones(1,length(constraints));
NewConstraints = [];
for i = 1:length(constraints)
    if is(constraints(i),'sdp')
        B0 = sdpvar(constraints(i));
        S = spy(B0);
        p = symrcm(S);
        bw = yalmipbandwidth(S(p,p));
        n = length(B0);
        if bw/n <= ratio & n>nmin
            kept(i) = 0;
            mid = floor(n/2+bw/2);
            B0 = B0(p,p);
            bwh = ceil(bw/2);
            NewConstraints = [NewConstraints, B0(1:mid,1:mid) >= 0, B0(mid-bw+1:end,mid-bw+1:end)];
        end
    end
end
if ~isempty(NewConstraints)
    NewConstraints = unblkdiag(NewConstraints);
    NewConstraints = splitAndDiagonalize(NewConstraints,nmin,ratio);
end
NewConstraints = [NewConstraints, constraints(find(kept))];
    