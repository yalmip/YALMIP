function p = addDisjointBilinearSDPcut(p)

if any(p.K.s) && any(p.variabletype == 1)
    top = 1+p.K.f+p.K.l+sum(p.K.q)+p.K.e+sum(p.K.p);
    newData = [];
    newKs = [];
    for i = 1:length(p.K.s)
        data = p.F_struc(top:top+p.K.s(i)^2-1,:);
        used = find(any(data(:,2:end)));
        if ~any(p.variabletype(used))
            % This is a linear SDP F(x) >= 0            
            possible_t = setdiff(p.linears,used);
            weCanDoIt = 1;
            for j = 1:length(possible_t)
                t = possible_t(j);                
                generateMonomials = t_times_x(t,used,p);
                if ~any(isnan(generateMonomials))
                    if ~isinf(p.lb(t))
                        % F0*(t-L) + Fi*xi*(t-L) >= 0
                        temp = data;
                        temp(:,1 + t) = temp(:,1);
                        temp(:,1) = temp(:,1)*-p.lb(t);
                        temp(:,1 + generateMonomials) = temp(:,1+used);
                        temp(:,1+used) = temp(:,1+used)*-p.lb(t);
                        newData = [newData;temp];
                        newKs = [newKs p.K.s(i)];
                    end
                    if ~isinf(p.ub(t))
                        % F0*(U-t) + Fi*xi*(U-t) >= 0
                        temp = data;
                        temp(:,1 + t) = -temp(:,1);
                        temp(:,1) = temp(:,1)*p.ub(t);
                        temp(:,1 + generateMonomials) = -temp(:,1+used);
                        temp(:,1+used) = temp(:,1+used)*p.ub(t);
                        newData = [newData;temp];
                        newKs = [newKs p.K.s(i)];
                    end
                end
            end
        end
        top = top + p.K.s(i)^2;
    end 
    if ~isempty(newData);
        p.F_struc = [p.F_struc;newData];
        p.K.s = [p.K.s newKs];
    end
end
   
function index = t_times_x(t,x,p)
index = [];
for i = 1:length(x)
    xi = x(i);    
    p_t = p.monomtable(t,:);
    p_x = p.monomtable(xi,:);    
    j = findrows(p.monomtable,p_t+p_x);
    if isempty(j)
        index(i) = nan;
    else
        index(i) = j;
    end
end
