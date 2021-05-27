function vexity = DeriveVexityFromInflection(properties,xL,xU)

% [point (changes to convex = 1/changes to concave = -1) ... ]
if isa(properties.inflection,'function_handle')
    data = properties.inflection(xL,xU);
    if isempty(data)
        vexity = 'none';
        return
    end
	inflections = [-inf data(1:2:end) inf];
else
    data = properties.inflection;
    inflections = [-inf data(1:2:end) inf];
end
convex = [-data(2:2:end)  data(end)];
vexity = 'none';
for k = 1:length(inflections)
    if xL >= inflections(k) && xU <= inflections(k+1)        
        if convex(k) == 1
            vexity = 'convex';
            return
        else
            vexity = 'concave';
            return
        end
    end
end