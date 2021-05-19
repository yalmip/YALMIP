function model = propagateInitial(model,tempF,tempK,tempx0)

if any(isnan(model.x0)) & ~all(isnan(model.x0))
    xevaled = zeros(1,length(model.c));
    xevaled(model.linearindicies) = model.x0;   
    xevaledout = apply_recursive_evaluation(model,xevaled);
    startNan = nnz(isnan(xevaledout));
    goon = 1;
    while goon
        temp = propagatex0(xevaledout,tempF,tempK);
        model.x0 = temp(model.linearindicies);
        if any(isnan(model.x0))
            xevaled(model.linearindicies) = model.x0;
            xevaledout = apply_recursive_evaluation(model,xevaled);
            goon =  nnz(isnan(xevaledout)) < startNan;
            startNan = nnz(isnan(xevaledout));
        else
            goon = 0;
        end
    end
end
model = correctEXPConeClosureInitial(model);
model.x0(isnan(model.x0))=0;