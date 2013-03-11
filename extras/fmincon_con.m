function [g,geq,dg,dgeq,xevaled] = fmincon_con(x,model,xevaled)

global latest_xevaled
global latest_x_xevaled
% Early bail for linear problems
g = [];
geq = [];
dg = [];
dgeq = [];
if model.linearconstraints
    xevaled = [];
    return
end

if nargin<3
    if isequal(x,latest_x_xevaled)
        xevaled = latest_xevaled;
    else
        xevaled = zeros(1,length(model.c));
        xevaled(model.linearindicies) = x;
        xevaled = apply_recursive_evaluation(model,xevaled);
        latest_x_xevaled = x;
        latest_xevaled = xevaled;
    end
end

if model.nonlinearinequalities
    g = full(model.Anonlinineq*xevaled(:)-model.bnonlinineq);
end

if model.nonlinearequalities
    geq = full(model.Anonlineq*xevaled(:)-model.bnonlineq);
end

dgAll_test = [];

if nargout == 2 || ~model.derivative_available
    return
elseif ~isempty(dgAll_test) & isempty(model.evalMap)
    dgAll = dgAll_test;
elseif isempty(model.evalMap) & (model.nonlinearinequalities | model.nonlinearequalities)
   % allA = [model.Anonlineq;model.Anonlinineq];
   % dgAll = [];
    n = length(model.c);
    linearindicies = model.linearindicies;
    %mtNonlinear = model.monomtable(model.nonlinearindicies,:);
    xevaled = zeros(1,n);
    xevaled(linearindicies) = x;
    % FIXME: This should be vectorized
    
    news = model.fastdiff.news;
    allDerivemt = model.fastdiff.allDerivemt;
    c = model.fastdiff.c;
    
    if model.fastdiff.univariateDifferentiates
        zzz = c.*(x(model.fastdiff.univariateDiffMonom).^model.fastdiff.univariateDiffPower);
    else
        %  X = repmat(x(:)',length(c),1);
        O = ones(length(c),length(x));
        nz = find(allDerivemt);
        %  O(nz) = X(nz).^allDerivemt(nz);
        O(nz) = x(ceil(nz/length(c))).^allDerivemt(nz);        
        zzz = c.*prod(O,2);               
    end
    
    if 1
        if 0
        a1 = news(1:length(c),2);
        a2 = model.nonlinearindicies(news(1:length(c),1))';
        nn = max(max(length(linearindicies)),max(news(:,2)));
        mm = max(max(linearindicies),max(model.nonlinearindicies(news(:,1))));
        %  newdxx = spalloc(length(linearindicies),max(a2),length(linearindicies));
        %   newdxx = spalloc(nn,mm,length(linearindicies));
        %   iii = sub2ind(size(newdxx),a1(:),a2(:));
        %   newdxx(iii) = zzz;
        
        % Moved from the for-loop below*
        a1 = [a1(:);(1:length(linearindicies))'];
        a2 = [a2(:);linearindicies(:)];
        %zzz = [zzz;repmat(1,length(linearindicies),1)];
        zzz = [zzz;ones(length(linearindicies),1)];
        newdxx = sparse(a1,a2,zzz,nn,mm);
        else
            newdxx = model.fastdiff.newdxx;
            newdxx(model.fastdiff.linear_in_newdxx) = zzz;
        end
        
        %    newdxx = spalloc(length(linearindicies),max(a2),length(linearindicies));
        %    iii = sub2ind(size(newdxx),a1,a2);
        %    newdxx(iii) = zzz;
        %    newdxx = sparse(a1,a2,zzz);
    else
        newdxx = spalloc(length(linearindicies),max(linearindicies),length(linearindicies));
        for i = 1:length(c)
            newdxx(news(i,2),model.nonlinearindicies(news(i,1)))=zzz(i);%c(i)*prod(x(:)'.^allDerivemt(i,:));
        end
    end

    % * moved up
    %ii = sub2ind(size(newdxx),1:length(linearindicies),linearindicies);
%     %newdxx(ii) = 1;
%     for i = 1:length(linearindicies)
%        newdxx(i,linearindicies(i)) = 1;
%     end
    
    newdxx = newdxx';
    if ~isempty(model.Anonlineq)
    %    newdxx = newdxx';
        dgAll = model.Anonlineq*newdxx;
    else
        dgAll = [];
    end
    if ~isempty(model.Anonlinineq)
    %    newdxx = newdxx';
        aux = model.Anonlinineq*newdxx;
        dgAll = [dgAll;aux];
    end    
  %  dgAll = allA*newdxx';
    
else
    allA = [model.Anonlineq;model.Anonlinineq];
  %  requested = any(allA',2);
  %  [i,j,k] = find((model.deppattern(find(requested),:)));
  %  requested(j) = 1;
    requested = model.fastdiff.requested;
    dx = apply_recursive_differentiation(model,xevaled,requested,model.Crecursivederivativeprecompute);
    dgAll = allA*dx;
end

if  model.nonlinearequalities
    dgeq = dgAll(1:size(model.Anonlineq,1),:)';
end
if model.nonlinearinequalities
    dg = dgAll(size(model.Anonlineq,1)+1:end,:)';
end
