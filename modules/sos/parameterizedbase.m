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


    if 1
        [ii,vars] = find(parametric_basis);
        ii = ii(:)';
        vars = vars(:)';
    else
        ii = [];
        vars = [];
        js = sum(parametric_basis,1);
        indicies = find(js);
        for i = indicies
            if js(i)
                j = find(parametric_basis(:,i));
                ii = [ii j(:)'];
                vars = [vars repmat(i,1,js(i))];
            end
        end
    end

    k = setdiff1D(1:n,ii);
    if isempty(k)
        p_base_parametric = p_base_parametric.*sparse(ii,repmat(1,1,n),params(vars));
    else
        pp = params(vars); % Must do this, bug in ML 6.1 (x=sparse(1);x([1 1]) gives different result in 6.1 and 7.0!)
        p_base_parametric = p_base_parametric.*sparse([ii k(:)'],repmat(1,1,n),[pp(:)' ones(1,1,length(k))]);
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
