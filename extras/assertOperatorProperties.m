function properties = assertOperatorProperties(properties)

properties = assertProperty(properties,'definiteness','none');
properties = assertProperty(properties,'convexity','none');
properties = assertProperty(properties,'monotonicity','none');
properties = assertProperty(properties,'smoothness',inf);
properties = assertProperty(properties,'symmetry','none');
properties = assertProperty(properties,'periodic',[]);
properties = assertProperty(properties,'shape','none');
properties = assertProperty(properties,'function',[]);
properties = assertProperty(properties,'derivative',[]);
properties = assertProperty(properties,'inverse',[]);
properties = assertProperty(properties,'convexhull',[]);
properties = assertProperty(properties,'bounds',[]);
properties = assertProperty(properties,'inversebounds',[]);
properties = assertProperty(properties,'f_upper',[]);
properties = assertProperty(properties,'f_lower',[]);
properties = assertProperty(properties,'df_upper',[]);
properties = assertProperty(properties,'df_lower',[]);
properties = assertProperty(properties,'domain',[-inf inf]);
properties = assertProperty(properties,'stationary',[]);
properties = assertProperty(properties,'inflection',[]);
properties = assertProperty(properties,'singularity',[]);
properties = assertProperty(properties,'discontinuity',[]);
properties = assertProperty(properties,'forbidden',[]);
properties = assertProperty(properties,'degree',[]);
properties = assertProperty(properties,'replace',[]);
if isa(properties.definiteness,'char')
    switch properties.definiteness
        case 'positive'
            properties = assertProperty(properties,'range',[0 inf]);
            properties.range(1) = max(properties.range(1),0);
        case 'negative'
            properties = assertProperty(properties,'range',[-inf 0]);
            properties.range(2) = min(properties.range(2),0);
        otherwise
            properties = assertProperty(properties,'range',[-inf inf]);
    end
else
    properties = assertProperty(properties,'range',[-inf inf]);
end
if properties.range(1) >= 0;properties.definiteness ='positive';end
if properties.range(2) <= 0;properties.definiteness ='negative';end

if isequal(properties.shape,'s-shape')
    properties.monotonicity = 'increasing';
end
if isequal(properties.shape,'z-shape')
    properties.monotonicity = 'decreasing';
end

properties = assertProperty(properties,'model','unspecified');

if isequal(properties.convexity,'nonconvex')
    properties.convex = 'none';
end
if isequal(properties.definiteness,'indefinite')
    properties.definiteness = 'none';
end
if isequal(properties.monotonicity,'nonmonotone')
    properties.monotonicity = 'none';
end
if isequal(properties.smoothness,'smooth')
    properties.smoothness = inf;
elseif isequal(properties.smoothness,'discontinuous')
    properties.smoothness = 0;
elseif isequal(properties.smoothness,'differentiable')
    properties.smoothness = 1;
end