function model = update_monomial_bounds(model,these)

if nargin == 1 & all(model.variabletype<=2) & any(model.variabletype)
    % Fast code for purely quadratic case
    x = model.bilinears(:,2);
    y = model.bilinears(:,3);
    z = model.bilinears(:,1);
    corners = [model.lb(x).*model.lb(y) model.ub(x).*model.lb(y) model.lb(x).*model.ub(y) model.ub(x).*model.ub(y)];
    
    maxz = max(corners,[],2);
    minz = min(corners,[],2);
    
    model.lb(z) = max(model.lb(z),minz);
    model.ub(z) = min(model.ub(z),maxz);
    return
end

if nargin == 1
    polynomials = find((model.variabletype ~= 0));
else
    polynomials = find((model.variabletype ~= 0));
    polynomials = polynomials(find(ismember(polynomials,these)));
end

for i = 1:length(polynomials)
    j = polynomials(i);
    if j<=length(model.lb)
        monomials = model.monomtable(j,:);
        bound = powerbound(model.lb,model.ub,monomials);
        model.lb(j) = max(model.lb(j),bound(1));
        model.ub(j) = min(model.ub(j),bound(2));
        [inversebound,var] = inversepowerbound(model.lb,model.ub,monomials, polynomials(i));
        if ~isempty(var)
            model.lb(var) = max(model.lb(var),inversebound(1));
            model.ub(var) = min(model.ub(var),inversebound(2));
        end
    end
end

function  [inversebound,var] = inversepowerbound(lb,ub,monomials,polynomial);
inversebound = [];
var = [];
[i,var,val] = find(monomials);
if all(val == fix(val)) & all(val >= 0)
    if length(var) == 1
        if even(val)
            if val > 2
                inversebound = [-inf inf];
                aux = inf;
                if ~isinf(lb(polynomial))
                    if lb(polynomial) >= 0
                        aux = lb(polynomial)^(1/val);
                    end
                end
                if ~isinf(ub(polynomial))
                    if ub(polynomial) >= 0
                        aux = max(aux,ub(polynomial)^(1/val));
                    end
                else
                    aux = inf;
                end
                inversebound = [-aux aux];
            else
                var = [];
            end
        elseif val >= 3
            inversebound = [-inf inf];
            if ~isinf(lb(polynomial))
                inversebound(1,1) = sign(lb(polynomial))*abs(lb(polynomial))^(1/val);
            end
            if ~isinf(ub(polynomial))
                inversebound(1,2) = sign(ub(polynomial))*abs(ub(polynomial))^(1/val);
            end
        end
    else
        var = [];
    end
else
    var = [];
end


