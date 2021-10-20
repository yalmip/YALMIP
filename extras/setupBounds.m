function LU = setupBounds(F,options,extendedvariables)
nv = yalmip('nvars');
L = repmat(-inf,nv,1);
U = repmat(inf,nv,1);
yalmip('setbounds',1:nv,L,U);
% This is a hack to avoid bound propagation when this function is
% called from optimizer.m
if isfield(options,'avoidequalitybounds')
    LU = getbounds(lmi(F),0);
else
    LU = getbounds(lmi(F),[],[L U]);
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

% propagate bilinear. This is needed when implies etc are used in a
% bilinearly parameterized optimzer object
[mainmonomtable,variabletype] = yalmip('monomtable');
bilin = find(variabletype == 1);
if ~isempty(bilin)
    monomtable = mainmonomtable(bilin,:);
    [i,j] = find(monomtable');
    i = reshape(i,2,[]);
    x = i(1,:)';
    y = i(2,:)';
    z = bilin(:);
    lb = LU(:,1);
    ub = LU(:,2);
    corners = [lb(x).*lb(y) ub(x).*lb(y) lb(x).*ub(y) ub(x).*ub(y)];
    
    maxz = max(corners,[],2);
    minz = min(corners,[],2);
    
    LU(bilin,1) = max(LU(bilin,1),minz);
    LU(bilin,2) = min(LU(bilin,2),maxz);
end

% propagate simple quadratic
quadratic = find(variabletype == 2);
if ~isempty(quadratic)
    monomtable = mainmonomtable(quadratic,:);    
    [i,j] = find(monomtable');   
    x = i;    
    z = quadratic(:);
    lb = LU(:,1);
    ub = LU(:,2);
    corners = [lb(x).*lb(x) ub(x).*lb(x) lb(x).*ub(x) ub(x).*ub(x)];
    
    maxz = max(corners,[],2);
    minz = max(min(corners,[],2),0);
    
    LU(quadratic,1) = max(LU(quadratic,1),minz);
    LU(quadratic,2) = min(LU(quadratic,2),maxz);
end
b = yalmip('tempbinvariables');
if ~isempty(b)
    LU(b,1) = max(0,LU(b,1));
    LU(b,2) = min(1,LU(b,2));
end
yalmip('setbounds',1:nv,LU(:,1),LU(:,2));