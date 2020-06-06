function vexity = DeriveVexityFromInflection(properties,xL,xU)

% [point (changes to convex = 1/changes to concave = -1) ... ]
inflections = [-inf properties.inflection(1:2:end) inf];
vexity = 'none';
for k = 1:length(inflections)-1
    if xL >= inflections(k) && xU <= properties.inflection(k+1)        
        if properties.inflection(2*k) == -1
            vexity = 'convex';
        else
            vexity = 'concave';
        end
    end
end