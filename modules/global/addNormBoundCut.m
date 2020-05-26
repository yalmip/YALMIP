function pcut = addNormBoundCut(p)

if isempty(p.bilinears) | p.K.l == 0 | ~any(p.variabletype == 2)
    pcut = p;
else
    pcut = emptyNumericalModel;
    
    % Search for (a_i >0)*x_i^2 >= bi, xi >= 0
    for i = p.K.f+1:p.K.f + p.K.l
        row = p.F_struc(i,:);
        b = row(1);
        if b < 0
            a = row(2:end);
            j = find(a);
            if all(a(j)>0) && all(p.variabletype(j)==2)
                % Found it!
                [~,loc] = ismember(j,p.bilinears(:,1));
                vars = p.bilinears(loc,2);
                if all(p.lb(vars)>=0)
                    % We have a1x1^2 + ... + anxn^2 >= bn
                    % i.e. (a1^0.5*x1+.. + an^0.5)^2 >= bn+sum cijxixj
                    % all coefficients and x positive, so it must be
                    % i.e. (a1^0.5*x1+.. + an^0.5) >= sqrt(bn+prod(lb))
                    acut = a*0;acut(vars) = sqrt(a(j));
                    pcut.F_struc = [pcut.F_struc; -sqrt(-b + prod(p.lb(vars))) acut];
                    pcut.K.l = pcut.K.l + 1;
                end
            end
        end
    end        
    pcut = mergeNumericalModels(p,pcut);
end



