function p_lp = addTrivialActivationCuts(p,p_lp)
if any(p.K.s)
    top = startofSDPCone(p.K);
    for k = 1:length(p.K.s)
        index = top:top+p.K.s(k)^2-1;
        F0 = p.F_struc(top:top+p.K.s(k)^2-1,1);
        for col = 1:p.K.s(k)
            for row = 2:p.K.s(k)
                pos = (col-1)*p.K.s(k) + row - 1 + top;
                if p.F_struc(pos,1) ~= 0
                    if ~any(p.F_struc(pos,2:end))
                        % Constant in (row,row)
                        pos = p.K.s(k)*(row-1)+row-1+top;
                        if p.F_struc(pos,1)==0
                            used = find(p.F_struc(pos,2:end));
                            if all(ismember(used,p.binary_variables))
                            elseif all(ismember(used,p.integer_variables))
                                if all(p.lb(used)>=0)
                                elseif all(p.ub(used)<=0)
                                    % All negative integers. Some of them 
                                    % has to be -1 at least
                                    p_lp.F_struc(end+1,1)=-1;
                                    p_lp.F_struc(end,used+1) = -1;
                                    p_lp.K.l = p_lp.K.l + 1;
                                end   
                            end
                        end
                    end
                end
            end                                       
        end
    end
end