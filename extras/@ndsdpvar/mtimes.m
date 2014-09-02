function F = mtimes(X,Y)
% mtimes (overloaded)

Xd = isa(X,'ndsdpvar');
Yd = isa(Y,'ndsdpvar');

if Xd & Yd
    error('nD SDPVAR objects can not be multiplied');
end

if ~Xd
    if isa(X,'double')
        if prod(size(X)) == 1
            F = Y;
            F.basis = F.basis*X;
            F = flush(F);
        else
            error('Only scalar multiplication allowed on nD objects');
        end
    end
end

if ~Yd
    if isa(Y,'double')
        if prod(size(Y)) == 1
            F = X;
            F.basis = F.basis*Y;
            F = flush(F);            
        else
            error('Only scalar multiplication allowed on nD objects');
        end
    end
end