function problematicQP = nonconvexQuadratic(Q)
problematicQP=0;
% Try to figure out if the quadratic part is convex or not.
% Avoid eigenvalue check etc as far as possible. This matrix might be huge!
[i,j,s] = find(Q);
if all(i == j)
    % diagonal
    if any(s<-1e-14)
        problematicQP = 1;
    end
    return
end

% Remove zero columns/rows
rmv = find(~any(Q,1));
Q(:,rmv)=[];
Q(rmv,:)=[];
[R,p] = chol(Q);
if ~p
    return
end

zeroDiag = find(diag(Q)==0);
if ~isempty(zeroDiag)
    if any(Q(:,zeroDiag))
        problematicQP = 1;
        return
    end
end

% A common case is that the matrix is block diagonal. Hence, we detect the
% blocks, and then check the blocks
[p,q,r,s,cc,rr] = dmperm(Q);

if ~all(p==q)
    problematicQP = 1;
    return
end

Q = Q(p,p);
block = 1;
while block <= length(r)-1
    
    if r(block)+1 == r(block+1)
        % Go thorough sequence of diagonal blocks
        diagend = block + 1;
        while diagend <= length(r)-1 & (r(diagend)+1 == r(diagend+1))
            diagend = diagend + 1;
        end
        diagend = diagend-1;
        Qblock = Q(r(block):r(diagend),r(block):r(diagend));
        if any(diag(Qblock)<-1e-14)
            problematicQP = 1;
            return
        end
        block = diagend+1;
    else
        Qblock = Q(r(block):r(block+1)-1,r(block):r(block+1)-1);
        
        block = block + 1;
        if length(Qblock)==1
            if Qblock < 1e-14
                problematicQP = 1;
                return
            end
        else
            [R,p] = chol(Qblock);
            if p
                if min(eig(full(Qblock)))<-1e-14
                    problematicQP = 1;
                    return
                end
            end
        end
    end
end