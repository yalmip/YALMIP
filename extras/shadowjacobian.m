function dfdx = shadowjacobian(f,x)
% See SDPVAR/jacobian

if isa(f,'double')
    dfdx = zeros(size(f,1),length(x));
    return
end

% if ~isempty(intersect(deepdepends(f),depends(x)))    
%     % Under development   
% end

if nargin==1
    if isa(f,'sdpvar')
        x = recover(depends(f));
    else
        x = 0;
    end
else
    if length(getvariables(x))<length(x)
      error('x should be a vector of scalar independant variables');
    end
end

[n,m]=size(f);

if m>1
   error('Jacobian only defined for column vectors.')
end

if n*m==1
    dfdx = scalar_jacobian(f,x);
    % Argh, fix this (sorts inside scalar_jacobian    
    for i = 1:length(x)
        var(i,1)=getvariables(x(i));
    end
    [i,j]=sort(var);
    dfdx = dfdx(1,j);
    return
    
else
    dfdx = [];
    AllVars = recover(unique([depends(f) getvariables(x)]));

    for i = 1:length(f)
        dfdx = [dfdx;scalar_jacobian(f(i),x,AllVars)];
    end
    % Argh, fix this (sorts inside scalar_jacobian
    for i = 1:length(x)
        var(i,1)=getvariables(x(i));
    end
    [i,j]=sort(var);
    dfdx = dfdx(:,j);
end

function [dfdx,dummy] = scalar_jacobian(f,x,AllVars)

if isa(f,'double')
    dfdx = zeros(1,length(x));
    return
end

if nargin==2
    AllVars = recover(uniquestripped([depends(f) getvariables(x)]));
    %AllVars = recover(uniquestripped([deepdepends(f) depends(f) getvariables(x)]));
end



exponent_p = exponents(f,AllVars);

coefficients = getbase(f);
coefficients = coefficients(2:end);
coefficients = coefficients(:);
if nnz(exponent_p(1,:))==0
    exponent_p=exponent_p(2:end,:);
end

x_variables = getvariables(x);
AllVars_variables = getvariables(AllVars);
%AllDeriv = [];
AllDeriv2 = [];
for k = 1:length(x)
    wrt = find(ismembcYALMIP(AllVars_variables,x_variables(k)));
    deriv = exponent_p;
    deriv(:,wrt) = deriv(:,wrt)-1;     
    keep{k} = find(deriv(:,wrt)~=-1);
    AllDeriv2 = [AllDeriv2;deriv(keep{k},:)];        
end

if size(AllDeriv2,1)==0
    dummy = 1;
else
    dummy = recovermonoms(AllDeriv2,AllVars);
end

top = 1;
dfdx=[];
for k = 1:length(x)
    wrt = find(ismembcYALMIP(AllVars_variables,x_variables(k)));    
    m = length(keep{k});
    if m>0
        poly = sum((coefficients(keep{k}(:)).*exponent_p(keep{k}(:),wrt)).*dummy(top:top+m-1),1);
        top = top + m;
    else
        poly = 0;
    end
    dfdx = [dfdx ; poly];
end
dfdx = dfdx';
