function P = CreateBasicOperator(varargin)

for i = 1:length(varargin)
    
    k = strmatch_octavesafe(varargin{i}, {'convex','concave','nonconvex'});
    if ~isempty(k)
        P.convexity = varargin{i};
        continue
    end
    
    k = strmatch_octavesafe(varargin{i}, {'increasing','decreasing','nonmonotone'});
    if ~isempty(k)
        P.monotonicity = varargin{i};
        continue
    end
    
    k = strmatch_octavesafe(varargin{i}, {'positive','negative','indefinite'});
    if ~isempty(k)
        P.definiteness = varargin{i};
        continue
    end
    
    k = strmatch_octavesafe(varargin{i}, {'exact','milp','callback','graph'});
    if ~isempty(k)
        P.model = varargin{i};
        continue
    end
    
end