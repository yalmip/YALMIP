function keep = newtonpolytope(exponent_m,exponent_p)
%NEWTONPOLYTOPE  Internal function to remove monimials in SOS programs using Newton polytope

% Author Johan Löfberg
% $Id: newtonpolytope.m,v 1.1 2006-03-30 13:27:20 joloef Exp $


% WARNING : THIS CODE SUCKS AND IS ONLY USED AS A BACK UP PLAN 
% IF EVERYTHING ELSE FAILS (CRASHING LP SOLVERS ETC)

% *************************************
% TRY TO CALCULATE CONVEX HULL
% *************************************
try
    cnvhull = convhulln(full(exponent_p));
catch
    keep = 1:size(exponent_m,1);
    return
end

% ***************************************
% GET THE UNIQUE POINTS OF IN CONVEX HULL
% ***************************************
unique_points = unique(cnvhull);

% ***************************************
% CALCULATE A POINT IN INTERIOR
% ***************************************
p_c = sum(exponent_p(unique_points,:),1)'/length(unique_points);

% ***************************************
% CALCULATE HYPER-PLANES Ai^T(x-bi)=0
% ***************************************
j = 1;
A = [];
for i = 1:size(cnvhull,1)
    X = exponent_p(cnvhull(i,:),:)';
    y = X(:,1);
    dX = X(:,2:end)-repmat(X(:,1),1,size(X,2)-1);
    Atemp = null(full(dX'));
    % Is this a full-dimensional facet
    if size(Atemp,2)==1
        direction = (p_c-y)'*Atemp;
        if direction > 0
           A{j}=-Atemp;
       else
           A{j}=Atemp;
       end
        b{j}=y;  
        j=j+1;
    end
end

% ***************************************
% CHECK IF MONOMIAL IS IN NEWTON POLYTOPE
% ***************************************
if ~isempty(A)
    keep = [];
    for j = 1:size(exponent_m,1)
        inside = 1;
        y = exponent_m(j,:)';
        if isempty(findrows(exponent_p,y)) % Numerically safe on border
            i = 1;
            while inside & (i<=length(A))
                inside = inside & ((A{i}'*(y-b{i}))<=1e-9);
                i = i+1;
            end
        end
        if inside
            keep = [keep j];
        end
    end
else
    keep = 1:size(exponent_m,1);
end
