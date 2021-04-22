function p = detectMonotoneSDPObjective(p)
p.monotoneobjectiveresponse = 0;
if isempty(p.evalMap) && nnz(p.Q) == 0 && all(p.variabletype == 0)
    s = setdiff(find(p.c),[p.binary_variables p.integer_variables]);
    if ~isempty(s)
        if all(p.c(s))>0
            if nnz(p.F_struc(1:p.K.f,1 + s))==0
                if all(all(p.F_struc(p.K.f+1:p.K.f+p.K.l,1 + s)>=0))
                    if all(all(p.F_struc(p.K.f+p.K.l+1:p.K.f+p.K.l+sum(p.K.q),1 + s)==0))
                        top = p.K.f + p.K.l + sum(p.K.q)+1;
                        working = 1;
                        for i = 1:length(p.K.s)
                            B = p.F_struc(top:top+p.K.s(i)^2-1,1+s);
                            % Variable enters at most one element (i.e. diagonal)
                            % and enters with positive sign
                            working = working==1 && all(sum(B | B,1)<=1) && all(all(B >=0));
                            top = top + p.K.s(i)^2;
                        end
                        if working
                            p.monotoneobjectiveresponse = 1;
                        end
                    end
                end
            end
        end
    end
end