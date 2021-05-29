function vexity = DeriveVexityFromInflection(properties,xL,xU)

% [point (changes to convex = 1/changes to concave = -1) ... ]
if isa(properties.inflection,'function_handle')
    data = properties.inflection(xL,xU);
    if isempty(data)
        vexity = 'none';
        return
    end
    data = [data inf];	
else
    data = properties.inflection;
    data = [data inf];    
end
convex = [data(2:2:end-1)];
points = data(1:2:end);
vexity = 'none';
for k = 1:length(points)-1
    if xL >= points(k) && xU <= points(k+1)        
        if convex(k) == 1
            vexity = 'convex';
            return
        else
            vexity = 'concave';
            return
        end
    end
end