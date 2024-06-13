function f = createConvexHullMethod(p,xL,xU)

if ~isa(p.properties.derivative,'function_handle')
    f = [];
    return
end

if isa(p.properties.convexity,'char') && ~(isequal(p.properties.convexity,'none'))
    vexity = p.properties.convexity;
elseif isa(p.properties.convexity,'function_handle')
    % User-supplied method to derive convexity in region
    vexity = p.properties.convexity(xL,xU);
elseif ~isempty(p.properties.inflection)
    % Derive convexity by information about inflection
    vexity = DeriveVexityFromInflection(p.properties,xL,xU);   
else
    vexity = 'none';
end

if isequal(vexity,'convex')
    if strcmpi(p.fcn,'blackbox')
        f0 = @(x)real(blackbox(x,p.arg{2}));
    else
        f0 = @(x)real(eval([p.fcn '(x)']));
    end
    f = @(xL,xU)createConvexHullMethodConvex(xL,xU,f0,p.properties.derivative);
elseif isequal(vexity,'concave')
    f0 = @(x)real(eval([p.fcn '(x)']));
    f = @(xL,xU)createConvexHullMethodConcave(xL,xU,f0,p.properties.derivative);
else
    f = [];
end

