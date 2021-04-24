function [base,v] = matrixcoefficients(p,x)
%MATRIXCOEFFICIENTS Extends coefficients, beta

% FIX: CURRENTLY UNSTABLE!

if nargout>1 & (max(size(p))>1)
    % error('For matrix inputs, only the coefficients can be returned. Request feature if you need this...');
end

% Hack to make sure we get the basis w.r.t all variables and 1
% This has to be fixed soon (to make robust opt. module fast)
p = p + pi + sum(x)*1e-5;

if nargin==1
    allvar = depends(p);
    xvar = allvar;
    x = recover(xvar);
else
    xvar = intersect(depends(x),depends(p));
end

% Try to debug this!
p = p(:);
base = [];
v = [];
allvar = depends(p);%(ii));
allvar_recovered = recover(allvar);
t = setdiff(allvar,xvar);
t_recovered = recover(t);
ParametricIndicies = find(ismember(allvar,t));
map = find(~ismember(allvar,t));

for ii = 1:length(p)
    pii = p(ii);
    [exponent_p,p_base] = getexponentbase(pii,allvar_recovered);

    tempbase = parameterizedbase(pii,[],t_recovered,ParametricIndicies,exponent_p,p_base);
    [i,j,k] = unique(full(exponent_p(:,map)),'rows');


    V = sparse(1:length(k),k,1,length(tempbase),max(k))';

    base{ii} = V*tempbase - [pi;repmat(1e-5,length(x),1)];

    keepthese = j(1:max(k));
    v{ii} = recovermonoms(exponent_p(keepthese,map),x);%recover(xvar));
end






function p_base_parametric = parameterizedbase(p,z, params,ParametricIndicies,exponent_p,p_base)

% Check for linear parameterization
parametric_basis = exponent_p(:,ParametricIndicies);
if all(sum(parametric_basis,2)==0)
    p_base_parametric = full(p_base(:));
    return
end
if all(sum(parametric_basis,2)<=1)
    p_base_parametric = full(p_base(:));
    n = length(p_base_parametric);
    ii = [];
    vars = [];
    js = sum(parametric_basis,1);
    for i = 1:size(parametric_basis,2)
        if js(i)
            j = find(parametric_basis(:,i));
            ii = [ii j(:)'];
            vars = [vars repmat(i,1,js(i))];
        end
    end
    k = setdiff1D(1:n,ii);
    if isempty(k)
        p_base_parametric = p_base_parametric.*sparse(ii,repmat(1,1,n),params(vars));
    else
        pp = params(vars); % Must do this, bug in ML 6.1 (x=sparse(1);x([1 1]) gives different result in 6.1 and 7.0!)
        temp = sparse([ii k(:)'],repmat(1,1,n),[pp(:)' ones(1,1,length(k))]);
        p_base_parametric = p_base_parametric.*temp;
    end
else
    % Bummer, nonlinear parameterization sucks...
    for i = 1:length(p_base)
        j = find(exponent_p(i,ParametricIndicies));
        if ~isempty(j)
            temp = p_base(i);
            for k = 1:length(j)
                if exponent_p(i,ParametricIndicies(j(k)))==1
                    temp = temp*params(j(k));
                else
                    temp = temp*params(j(k))^exponent_p(i,ParametricIndicies(j(k)));
                end
            end
            xx{i} = temp;
        else
            xx{i} = p_base(i);
        end
    end
    p_base_parametric = stackcell(sdpvar(1,1),xx)';
end
