function p = presolve_implied_binaryproduct(p)
% presolve y = x1*x2 where x1 x2 binary
% thus y is binary (and y = x1&x2 which might be exploited)
% Modelled by y >= x1 + x2 - 1, x1>=y, x2>=y
p.binaryProduct = [];
if p.K.l == 0
    return
end
% Extract linear inequalities
ff = p.F_struc(startofLPCone(p.K):startofLPCone(p.K)+p.K.l-1,:);
% Search for y >= x1 + x2 - 1
candidates = find(sum(abs(ff),2) == 4 & sum(ff,2) == 0 & ff(:,1)==1);
lower_and = {};
for i = 1:length(candidates)
    row = p.F_struc(candidates(i),2:end);
    xx = find(row==-1);
    yy = find(row==1);    
    if length(xx)==2 && length(yy)==1
        % Yes, this matches y-x1-x2+1>=0
        if all(ismember(xx,p.binary_variables))
            lower_and{end+1}.x = xx;
            lower_and{end}.y = yy;
        end
    end
end
for i = 1:length(lower_and)
    x = lower_and{i}.x;
    y = lower_and{i}.y;
    y_bounded_by_x1 = 0;
    y_bounded_by_x2 = 0;
    % Generalization of y <= x1 in some models
    candidates = find(ff(:,1) == 0 & ff(:,x(1)+1)==1 & ff(:,y+1)==-1);
    for j = 1:length(candidates)
        row = ff(candidates(j),2:end);
        row(x(1))=0;        
        % y + sum z <= x1   (z non-negative)
        % i.e x1 - sum(z) - y >= 0
        if all(row<=0) && all(p.lb(find(row))>=0)
            y_bounded_by_x1 = 1;
            break
        end
    end
    if y_bounded_by_x1
        candidates = find(ff(:,1) == 0 & ff(:,x(2)+1)==1 & ff(:,y+1)==-1);
        for j = 1:length(candidates)
            row = ff(candidates(j),2:end);
            row(x(2))=0;
            if all(row<=0) && all(p.lb(find(row)) >=0)
                y_bounded_by_x2 = 1;
                break
            end
        end
        if y_bounded_by_x2
            p.binary_variables = union(p.binary_variables,y);            
            p.binaryProduct = [p.binaryProduct;y x(:)']; 
            p.lb(p.binary_variables) = max(p.lb(p.binary_variables),0);
            p.ub(p.binary_variables) = min(p.ub(p.binary_variables),1);
        end
    end
end
p.integer_variables = setdiff(p.integer_variables,p.binary_variables);