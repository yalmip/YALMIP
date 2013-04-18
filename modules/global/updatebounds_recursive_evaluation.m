function p = updatebounds_recursive_evaluation(p)

if p.changedbounds
    if isempty(p.evalMap) & all(p.variabletype <= 2)
        % Bilinear/quadratic case can be done much faster
        p = updatemonomialbounds(p);
    else
        for i = 1:length(p.evaluation_scheme)
            switch p.evaluation_scheme{i}.group
                case 'eval'
                    for j = 1:length(p.evaluation_scheme{i}.variables)
                        p = update_one_eval_bound(p,j);
                        p = update_one_inverseeval_bound(p,j);
                    end
                case 'monom'
                    for j = 1:length(p.evaluation_scheme{i}.variables)
                        p = update_one_monomial_bound(p,j);
                    end
                otherwise
            end
        end
    end
    % This flag is turned on if a bound tightening funtion manages to
    % tighten the bounds
    p.changedbounds = 0;
end

function p = update_one_monomial_bound(p,indicies);
j = p.monomials(indicies);
bound = powerbound(p.lb,p.ub,p.monomtable(j,:));
p.lb(j) = max(p.lb(j),bound(1));
p.ub(j) = min(p.ub(j),bound(2));
