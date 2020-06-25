function monotonicity = DeriveMonotonicityFromShape(properties,xL,xU)

monotonicity = 'none';
if isequal(properties.shape,'bell-shape')
    if xU <= properties.stationary(1)
        monotonicity = 'increasing';
    elseif xL >= properties.stationary(1)
        monotonicity = 'decreasing';
    end
elseif isequal(properties.convexity,'v-shape')
    if xU <= properties.stationary(1)
        monotonicity = 'decreasing';
    elseif xL >= properties.stationary(1)
        monotonicity = 'increasing';
    end
end