function [M,m,infbound,lowerinf,upperinf] = derivebounds(f,lowerupper)
% Code to compute  estimates of max and min
% of the linear vector f
basis = getbase(f);
fvars = getvariables(f);
if nargin == 1
    lowerupper = yalmip('getbounds',fvars);
else
    lowerupper = lowerupper(fvars,:);
end

if any(any(isinf(lowerupper)))
    infbound = 1;
    badlower = isinf(lowerupper(:,1));
    lowerinf = any(badlower); 
    badupper = isinf(lowerupper(:,2));
    upperinf = any(badupper); 
    
    % Use artificial bound for unbounded variables
    % Make sure to be huge if model is badly scaled,
    % otherwise defult is 10^4
    largestU = max(lowerupper(~isinf(lowerupper(:))));
    smallestL = min(lowerupper(~isinf(lowerupper(:))));
    synthetic_bound = max(abs(largestU),abs(smallestL));
    synthetic_bound = max(10^4,synthetic_bound);
    if lowerinf && length(synthetic_bound)>0
        lowerupper(badlower,1) = -synthetic_bound;
    end
    if upperinf && length(synthetic_bound)>0
        lowerupper(badupper,2) = synthetic_bound;
    end

else
    infbound = 0;
    lowerinf = 0;
    upperinf = 0;
end

b = basis(:,1);
A = basis(:,2:end);
M = repmat(-inf,size(b,1),1);
m = -M;
if all(all(isinf(lowerupper)))
    M = repmat(inf,length(b),1);
    m = repmat(-inf,length(b),1);
else
    lower = lowerupper(:,1);
    upper = lowerupper(:,2);
    Apos = (A>0);
    Aneg = (A<0);
    npA = nnz(Apos); 
    nnA = nnz(Aneg);
    if all(all(~isinf(lowerupper)))
        M = A.*Apos*upper+A.*Aneg*lower;
        m = A.*Aneg*upper+A.*Apos*lower;
    elseif all(isinf(upper)) & nnA==0 & ~any(isinf(lower))
        % Speeds up norm(x,1) case
        m = A.*Apos*lower;    
    else
        for i = 1:length(b)
            ind = A(i,:)>0;
            i1 = find(ind);
            i2 = find(~ind);
            a = A(i,[i1 i2]);
            M(i) = a*[upper(i1);lower(i2)];
            m(i) = a*[lower(i1);upper(i2)];
        end
    end
end
M(isnan(M)) = inf;
M(isinf(M)) = 1e4; M = M+b;
m(isnan(m)) = -inf;
m(isinf(m)) = -1e4;m = m+b;
