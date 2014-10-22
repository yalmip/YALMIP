function NewConstraints = splitAndDiagonalize(constraints,nmin,ratio)

kept = ones(1,length(constraints));
NewConstraints = [];
for i = 1:length(constraints)
    if is(constraints(i),'sdp')
        B0 = sdpvar(constraints(i));
        S = spy(B0);
        p = symrcm(S);
        S = S(p,p);
        bw = yalmipbandwidth(S);
        n = length(B0);
        if bw/n <= ratio & n>nmin
            kept(i) = 0;
            mid = floor(n/2+bw/2);
            B0 = B0(p,p);
            bwh = ceil(bw/2);         
            E = sdpvar(bw+3).*S(mid-bw-1:mid+1,mid-bw-1:mid+1);
            B1 = B0(1:mid+1,1:mid+1);
            B1(end-bw-2:end,end-bw-2:end) = E;            
            B2 = B0(mid-bw-1:end,mid-bw-1:end);
            B2(1:bw+3,1:bw+3) =  B2(1:bw+3,1:bw+3)-E;
            NewConstraints = [NewConstraints, B1 >= 0,
                                              B2 >=0];
                                          
            %NewConstraints = [NewConstraints, B0(1:mid+1,1:mid+1)+blkdiag(zeros(2),-E) >= 0,
            %                                  B0(mid-bw-1:end,mid-bw-1:end)+blkdiag(E,zeros(3))>=0];
%             clf
%             H = zeros(n);
%             H(1:mid+1,1:mid+1)=1;
%             H(mid-bw-1:end,mid-bw-1:end)=1;
%             spy(H)
%             hold on;
%             spy(S(p,p),'r')
%             1;
        end
    end
end
if ~isempty(NewConstraints)
    NewConstraints = unblkdiag(NewConstraints);
    NewConstraints = splitAndDiagonalize(NewConstraints,nmin,ratio);
end
NewConstraints = [NewConstraints, constraints(find(kept))];
    