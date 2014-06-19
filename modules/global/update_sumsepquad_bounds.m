function p = update_sumsepquad_bounds(p);
% Looking for case z = b+ q1(x1) + q2(x2) + ... where q quadratic

if p.boundpropagation.sepquad
    found = 0;
    for j = 1:p.K.f
        a = p.F_struc(j,2:end);
        b = p.F_struc(j,1);
        used = find(a);
        data = [];
        if nnz(a) > 2 && all(p.variabletype(used) == 2 | p.variabletype(used) == 0)
            % Only linear or quadratic and there is a mixed term
            nonlinears = find(p.variabletype(used)==2);
            if (nnz(a) > length(nonlinears)) && length(nonlinears) > 0
                data = [];
                for i = 1:length(nonlinears)
                    linear = find(p.monomtable(nonlinears(i),:));
                    data = [data;linear a(linear) a(nonlinears(i))];
                    a(linear)=0;
                    a(nonlinears(i))=0;
                    found = 1;
                end
            end
        end
        if nnz(a) == 1 & ~isempty(data)
            k = find(a);
            ai = a(k);
            data(:,2:end) =  data(:,2:end)/(-ai);
            b = b/-ai;
            L = b;
            U = b;
            for i = 1:size(data,1)
                [Li,Ui] = wcquad(data(i,2:3),p.lb(data(i,1)),p.ub(data(i,1)));
                L = L + Li;
                U = U + Ui;
            end
            p.lb(k) = max(p.lb(k),L);
            p.ub(k) = min(p.ub(k),U);
        end
    end
    if ~found
        % Turn off this feature
        p.boundpropagation.sepquad = 0;
    end
end


function [Li,Ui] = wcquad(c,lb,ub)        

% Stationary point
xs = -c(1)/(2*c(2));
fs = c(1)*xs+c(2)*xs^2;
stationaryInBounds = (xs >= lb) && (xs <= ub);

if isinf(lb)
    if c(2)>0
        fl = inf;
    else
        fl = -inf;
    end
else
    fl = c(1)*lb+c(2)*lb^2;
end
if isinf(ub)
    if c(2)>0
        fu = inf;
    else
        fu = -inf;
    end
else
    fu = c(1)*ub+c(2)*ub^2;
end

if stationaryInBounds
    Li = min([fl fu fs]);
    Ui = max([fl fu fs]);
else
    Li = min([fl fu]);
    Ui = max([fl fu]);
end

if any(isnan([Li Ui]))
    error('hej')
end
         
   
   
 