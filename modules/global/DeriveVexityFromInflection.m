function vexity = DeriveVexityFromInflection(properties,xL,xU)
if xU <= properties.inflection(1)
    % Left region
    if properties.inflection(2) == -1
        vexity = 'convex';
    else
        vexity = 'concave';
    end
elseif xL >= properties.inflection(1)
    % Right region
    if properties.inflection(2) == -1
        vexity = 'concave';
    else
        vexity = 'convex';
    end
else
    % Covering inflection
    vexity = 'none';
end