function F = detset(t,P)
%DETSET Internal function used in construction of MAXDET formulations
%
% F = detset(t,P) creates the model t <= det(P)^(1/(2^ceil(log2(length(P)))))

[n,m]=size(P);

if max(n,m)==1
    F = (P>=0) + (t<=P);
else
    if min(n,m)==1
        % Vector version (copy and pasted from below)
        p = 2^ceil(log2(max(n,m)));
        x = [P(:);ones(p-max(n,m),1)];
        F = ([]);
    else
        % This code should never run!! (Taken care of in geomean.m)
        % Is P square?
        if n~=m
            error('P has to be square in the constraint t<det(P)^1/n')
        end
        % Is P symmetric?
        if ~is(P,'hermitian')
            error('P has to be Hermitian in the constraint t<det(P)^1/n')
        end
        % Is P complex?
        if is(P,'complex')
            P = [real(P) imag(P);-imag(P) real(P)];
            [n,m]=size(P);
        end

        D = tril(sdpvar(n,n));
        delta = diag(D);
        F = ([P D;D' diag(delta)] >= 0);
        p = 2^ceil(log2(n));
        x = [delta;t*ones(p-n,1)];
    end
    
    if 0
        % Slightly more efficient approach
        % Thanks to M Grant for pointing out this
        % (not thouroghly bug-tested with geomean yet)
        p = n;
        while p > 2
            p2 = ceil( p / 2 );
            if p2 * 2 ~= p,
                x = [ x ; t ];
            end
            x_new = sdpvar(p2,1);
            for i = 1:p2
                F = F + gmset(x_new(i),x(2*i-1),x(2*i));
            end
            x = x_new;
            p = p2;
        end
    else
        while p > 2
            x_new = sdpvar(p/2,1);
            for i = 1:p/2
                F = F + gmset(x_new(i),x(2*i-1),x(2*i));
            end
            x = x_new;
            p = p/2;
        end
    end
     F = F+gmset(t,x(1),x(2));
end
