function vexity = DeriveVexityFromInflection(properties,xL,xU)

% [point (changes to convex = 1/changes to concave = -1) ... ]
inflections = [-inf properties.inflection(1:2:end) inf];
convex = [-properties.inflection(2:2:end)  properties.inflection(end)];
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