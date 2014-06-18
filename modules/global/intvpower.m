function C = intvpower(A,p)
if p>0
    if nnz(rem(full(p),2))==0
        if A(1)<0 & A(2)>0
            C = [0 max([A(1)^p A(2)^p])];
        else
            a = A(1)^p;b = A(2)^p;
            C = [min([a b]) max([a b])];
        end
    else
        if  isequal(-inf,A(1)) && isequal(inf,A(2))
            C = A;
        elseif p~=fix(p) & A(1)<0
            C = [-inf inf];
        else
            a = A(1)^p;b = A(2)^p;
            C = [a b];
        end
    end
else
    if A(1)<0 && A(2)>0
        % Nasty crossing
        if even(p)
            C = [min([A(1)^p A(2)^p]) inf];
        else
            C = [-inf inf];
        end
    elseif A(1)>= 0
        % Decaying function
        C = [A(2)^p (abs(A(1)))^p];
    elseif even(p)
        
        a = A(1)^p;b = A(2)^p;
        C = [min([a b]) max([a b])];
    else
        a = A(1)^p;b = A(2)^p;
        C = [b a];
    end
end
