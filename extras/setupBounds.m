function LU = setupBounds(F,options,extendedvariables)
nv = yalmip('nvars');
yalmip('setbounds',1:nv,repmat(-inf,nv,1),repmat(inf,nv,1));
% This is a hack to avoid bound propagation when this function is
% called from optimizer.m
if isfield(options,'avoidequalitybounds')
    LU = getbounds(F,0);
else
    LU = getbounds(F);
end
% In models with nonconvex terms including x, but where bounds only are
% set on abs(x), we have to use the bound on abs(x) to improve the
% bound on x. Sort of ugly...The problem is that we only get operator
% knowledge when we model the operator, nor before we start the whole
% algorithm
% This should be generalized, operator.derivebounds
extstruct = yalmip('extstruct');
for i = 1:length(extstruct)
    switch extstruct(i).fcn
        case 'milpsubsref'
            LU = extract_bounds_from_milpsubsref_operator(LU,extstruct,extendedvariables,i);
        case 'abs'
            LU = extract_bounds_from_abs_operator(LU,extstruct,extendedvariables,i);
        case 'norm'
            LU = extract_bounds_from_norm_operator(LU,extstruct,extendedvariables,i);
        case 'min_internal'
            LU = extract_bounds_from_min_operator(LU,extstruct,extendedvariables,i);
        case 'max_internal'
            LU = extract_bounds_from_max_operator(LU,extstruct,extendedvariables,i);
        otherwise
    end
end
yalmip('setbounds',1:nv,LU(:,1),LU(:,2));