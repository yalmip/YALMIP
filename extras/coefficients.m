function [base,v] = coefficients(p,x,vin)
%COEFFICIENTS Extract coefficients and monomials from polynomials
%
%   [c,v] = COEFFICIENTS(p,x) extracts the coefficents
%   of a scalar polynomial p(x) = c'*v(x)
%
%   c = COEFFICIENTS(p,x) extracts the all coefficents
%   of a matrix polynomial.

%
%   INPUT
%    p : SDPVAR object
%    x : SDPVAR object
%
%   OUTPUT
%    c : SDPVAR object
%    v : SDPVAR object
%
%   EXAMPLE
%    sdpvar x y s t                
%    p = x^2+x*y*(s+t)+s^2+t^2;     % define p(x,y), parameterized with s and t
%    [c,v] = coefficients(p,[x y]); 
%    sdisplay([c v]) 
%
%   See also SDPVAR

% Author Johan Löfberg
% $Id: coefficients.m,v 1.19 2010-03-14 12:57:16 joloef Exp $  

if isa(p,'double')
    base = p(:);
    v = 1;
    return
end

if nargout>1 & (max(size(p))>1)
    error('For matrix inputs, only the coefficients can be returned. Request feature if you need this...');
end

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
for i = 1:length(p)
    allvar = depends(p(i));
    t = setdiff(allvar,xvar);
    if 1%isa(p(i),'sdpvar')
        [exponent_p,p_base] = getexponentbase(p(i),recover(depends(p(i))));
        ParametricIndicies = find(ismember(allvar,t));
        % FIX : don't define it here, wait until sparser below. Speed!!
        tempbase = parameterizedbase(p(i),[],recover(t),ParametricIndicies,exponent_p,p_base);
        [i,j,k] = unique(full(exponent_p(:,find(~ismember(allvar,t)))),'rows');
        V = sparse(1:length(k),k,1,length(tempbase),max(k))';
        base = [base;V*tempbase];
        if nargout == 2
            keepthese = j(1:max(k));
            v = recovermonoms(exponent_p(keepthese,find(~ismember(allvar,t))),recover(xvar));
        end
    end
    if isa(p,'ncvar')
        [exponent_p,ordered_list] = exponents(p,recover(depends(p(i))));
        ParametricIndicies = find(ismember(allvar,t));
        
        for j = 1:length(ParametricIndicies)
            a = find(ordered_list(:,1) == ParametricIndicies(j))
            b = [];
            for k = 1:length(a)
            b = [b ordered_list(a(k),2:end)]
            end
            b = b(find(b));
            basetemp = [];
            for k = 1:length(b)
                basetemp = [basetemp ncvar(struct(recover(t((k)))))];
            end
            base = [base;sum(basetemp)];
        end
    end
end

if nargout <= 1
    v = [];
    vin=v;
else
    if nargin<3
        vin=v;
    end
end
if isequal(v,vin)
    return
else
    for i = 1:length(v)
        if isa(v(i),'double')
            si(i) = 0;
        else
            si(i) = getvariables(v(i));
        end
    end
    for i = 1:length(vin)
        if isa(vin(i),'double')
            vi(i) = 0;
        else
            vi(i) = getvariables(vin(i));
        end
    end

    newcvals = [];
    if  all(ismember(si,vi))
        for i = 1:length(vin)
            where = find(vi(i) == si);
            if isempty(where)
                newcvals = [newcvals;0];
                %newc(i,1) = 0;
            else
                %newc(i,1) = base(where);
                newcvals = [newcvals;base(where)];
            end
        end
        newc = sparse(1:length(vin),ones(length(vin),1),newcvals);
    else
        error('The supplied basis is not sufficient');
    end

    base = newc(:);
    v = vin(:);
end


function p_base_parametric = parameterizedbase(p,z, params,ParametricIndicies,exponent_p,p_base)

% Check for linear parameterization
parametric_basis = exponent_p(:,ParametricIndicies);
%if all(sum(parametric_basis,2)==0)
if all(all(parametric_basis==0))
    p_base_parametric = full(p_base(:));
    return
end
if all(ismember(parametric_basis,[0 1])) & all(sum(parametric_basis,2)<=1)%all(sum(parametric_basis,2)<=1)
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



% 
% 
% 
% 
% 
% function [base,v] = coefficients(p,x)
% %COEFFICIENTS Extract coefficients and monomials from polynomials
% %
% %   [c,v] = COEFFICIENTS(p,x) extracts the coefficents
% %   of a polynomial p(x) = c'*v(x)
% %
% %   INPUT
% %    p : SDPVAR object
% %    x : SDPVAR object
% %
% %   OUTPUT
% %    c : SDPVAR object
% %    v : SDPVAR object
% %
% %   EXAMPLE
% %    sdpvar x y s t                
% %    p = x^2+x*y*(s+t)+s^2+t^2;     % define p(x,y), parameterized with s and t
% %    [c,v] = coefficients(p,[x y]); 
% %    sdisplay([c v]) 
% %
% %   See also SDPVAR
% 
% % Author Johan Löfberg
% % $Id: coefficients.m,v 1.19 2010-03-14 12:57:16 joloef Exp $  
% 
% 
% if length(p) > 1%size(p,2) > 1
%     error('Coefficents can only be applied to column vectors');
% end
% 
% allvar = depends(p);
% if nargin==1
%     xvar = allvar;
%     x = recover(xvar);
% else
%     xvar = depends(x);    
% end
% 
% pvar = recover(depends(p));
% 
% base = [];
% for i = 1:length(p)
%     [bi{i},vi{i}] = coefficientsi(p(i),xvar,pvar,allvar);
% end
% 
% % Fix the lengths of the basis to use same basis for all elements
% if length(bi)>1
%     allvars = [];
%     for i = 1:length(bi)
%         bivar{i} = getvariables(vi{i});
%         if isequal(vi{i}(1),1)
%             bivar{i} = [0 bivar{i}];
%         end
%         allvars = unique([allvars bivar{i}]);
%     end
%     v = recover(allvars);
%     c = zeros(length(p),length(allvars))';
%     ci = [];
%     cj = [];
%     cv = [];
%     for i = 1:length(bi)
%         index = find(ismember(allvars,bivar{i}));
%         ci = [ci index];
%         cj = [cj repmat(i,1,length(index))];
%         cv = [cv bi{i}'];        
%     end
%     base = sparse(ci,cj,cv);
% else
%     base = bi{1};
%     v = vi{1};
% end
% 
% 
% function [base,v] = coefficientsi(p,xvar,pvar,allvar)
% 
% % Try to debug this!
% t = setdiff(allvar,xvar);
% [exponent_p,p_base] = getexponentbase(p,pvar);
% ParametricIndicies = find(ismember(allvar,t));
% % FIX : don't define it here, wait until sparser below. Speed!!
% tempbase = parameterizedbase(p,[],recover(t),ParametricIndicies,exponent_p,p_base);
% [i,j,k] = unique(full(exponent_p(:,find(~ismember(allvar,t)))),'rows');
% %V = sparse(max(k),length(tempbase));
% %for i = 1:max(k)    
% %    V(i,find(k==i)) = 1;
% %end
% V = sparse(1:length(k),k,1,length(tempbase),max(k))';
% base = V*tempbase;
% if nargout == 2
%     keepthese = j(1:max(k));
%     v = recovermonoms(exponent_p(keepthese,find(~ismember(allvar,t))),recover(xvar));
% end
% 
% 
% function p_base_parametric = parameterizedbase(p,z, params,ParametricIndicies,exponent_p,p_base)
% 
% % Check for linear parameterization
% parametric_basis = exponent_p(:,ParametricIndicies);
% if all(sum(parametric_basis,2)==0)
%     p_base_parametric = full(p_base(:));
%     return
% end
% if all(sum(parametric_basis,2)<=1)
%     p_base_parametric = full(p_base(:));
%     n = length(p_base_parametric);
%     ii = [];
%     vars = [];
%     js = sum(parametric_basis,1);
%     for i = 1:size(parametric_basis,2)
%         if js(i)
%             j = find(parametric_basis(:,i));
%             ii = [ii j(:)'];
%             vars = [vars repmat(i,1,js(i))];
%         end
%     end
%     k = setdiff1D(1:n,ii);
%     if isempty(k)
%         p_base_parametric = p_base_parametric.*sparse(ii,repmat(1,1,n),params(vars));
%     else
%         pp = params(vars); % Must do this, bug in ML 6.1 (x=sparse(1);x([1 1]) gives different result in 6.1 and 7.0!)
%         p_base_parametric = p_base_parametric.*sparse([ii k(:)'],repmat(1,1,n),[pp(:)' ones(1,1,length(k))]);
%     end
% else
%     % Bummer, nonlinear parameterization sucks...
%     for i = 1:length(p_base)
%         j = find(exponent_p(i,ParametricIndicies));
%         if ~isempty(j)
%             temp = p_base(i);
%             for k = 1:length(j)
%                 if exponent_p(i,ParametricIndicies(j(k)))==1
%                     temp = temp*params(j(k));
%                 else
%                     temp = temp*params(j(k))^exponent_p(i,ParametricIndicies(j(k)));
%                 end
%             end
%             xx{i} = temp;
%         else
%             xx{i} = p_base(i);
%         end
%     end
%     p_base_parametric = stackcell(sdpvar(1,1),xx)';
% end
