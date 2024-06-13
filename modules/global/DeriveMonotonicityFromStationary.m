function monotonicity = DeriveMonotonicityFromStationary(properties,xL,xU)

monotonicity = 'none';
if isequal(properties.convexity,'convex')
    if xU <= properties.stationary(1)
        monotonicity = 'decreasing';
    elseif xL >= properties.stationary(1)
        monotonicity = 'increasing';
    end
elseif isequal(properties.convexity,'concave')
    if xU <= properties.stationary(1)
        monotonicity = 'increasing';
    elseif xL >= properties.stationary(1)
        monotonicity = 'decreasing';
    end
end