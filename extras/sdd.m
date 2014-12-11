function Constraint = sdd(X)

if issymmetric(X)
    Constraint = [];
    n = size(X,1);
    M = 0;
    for ii = 1:n
        for jj = [1:1:ii-1 ii+1:1:n]
           Mij = sdpvar(2);
           Constraint = Constraint + sdp2socp(Mij);
            M = M + sparse([ii jj ii jj],[ii ii jj jj],Mij(:),n,n);
        end
    end
    Constraint = Constraint + [M == X];
else
    error('sdd requires a symmetric argument.');
end

function F = sdp2socp(M)
F=rcone(M(1,2),.5*M(1,1),M(2,2));