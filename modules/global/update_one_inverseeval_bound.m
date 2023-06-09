function p = update_one_inverseeval_bound(p,i)

arg = p.evalMap{i}.variableIndex;
xL = p.lb(arg);
xU = p.ub(arg);
if ~all(xL==xU)
    fL = p.lb(p.evalMap{i}.computes);
    fU = p.ub(p.evalMap{i}.computes);
    try
        if ~isempty(p.evalMap{i}.properties.inverse)
            if strcmpi(p.evalMap{i}.properties.monotonicity,'increasing')
                invfiL = real(feval(p.evalMap{i}.properties.inverse,fL,p.evalMap{i}.arg{2:end-1}));
                invfiU = real(feval(p.evalMap{i}.properties.inverse,fU,p.evalMap{i}.arg{2:end-1}));
                p.lb(arg) = max(p.lb(arg),invfiL);
                p.ub(arg) = min(p.ub(arg),invfiU);
            elseif strcmpi(p.evalMap{i}.properties.monotonicity,'decreasing')
                invfiL = real(feval(p.evalMap{i}.properties.inverse,fU,p.evalMap{i}.arg{2:end-1}));
                invfiU = real(feval(p.evalMap{i}.properties.inverse,fL,p.evalMap{i}.arg{2:end-1}));
                p.lb(arg) = max(p.lb(arg),invfiL);
                p.ub(arg) = min(p.ub(arg),invfiU);
            end
        elseif ~isempty(p.evalMap{i}.properties.inversebounds)
            [invfiL,invfiU] = feval(p.evalMap{i}.properties.inversebounds,fL,fU,xL,xU);
             p.lb(arg) = max(p.lb(arg),invfiL);
             p.ub(arg) = min(p.ub(arg),invfiU);
        end
    catch
    end
end
