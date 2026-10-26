function n = sdpLayerConstraintCount(model)
% Number of scalar constraints used to represent every SDP block.

singleConstraintAlgorithm = ischar(model.options.slayer.algorithm) && ...
    any(strcmpi(model.options.slayer.algorithm,{'projector','signed-distance'}));
if singleConstraintAlgorithm
    n = ones(size(model.K.s));
elseif ischar(model.options.slayer.algorithm) && strcmpi(model.options.slayer.algorithm,'signed-distance-multi')
    n = repmat(length(model.options.slayer.p),size(model.K.s));
elseif isequal(model.options.slayer.m,inf)
    n = model.K.s;
elseif isequal(model.options.slayer.m,-1)
    n = min(length(model.c),model.K.s);
else
    n = min(model.K.s,repmat(model.options.slayer.m,1,length(model.K.s)));
end
