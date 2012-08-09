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
LU = extract_bounds_from_abs_operator(LU,yalmip('extstruct'),extendedvariables);
yalmip('setbounds',1:nv,LU(:,1),LU(:,2));