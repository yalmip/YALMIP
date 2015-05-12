function LU = getbounds(F,avoidequalitybounds,LU)

K.f = 0;
K.l = 0;
L = [];
U = [];
if nargin < 3
    LU =  yalmip('getbounds',1:yalmip('nvars'));
end
binary = yalmip('binvariables');
LU(binary,1) = 0;
LU(binary,2) = 1;
F = flatten(F);
simplex = [];
is_interval = is(F,'interval');
for i = 1:length(F.LMIid)
    if F.clauses{i}.type == 2
        X = F.clauses{i}.data;
        AB = getbase(X);
        K.l = prod(size(X));
        variables = getvariables(X);
        if is_interval(i)
            [lb,ub,cand_rows] = findulb_interval(AB,K);
        else
            [lb,ub,cand_rows] = findulb(AB,K);            
        end
        LU(variables,1) = max([lb LU(variables,1)]')';
        LU(variables,2) = min([ub LU(variables,2)]')';  
    elseif  F.clauses{i}.type == 55
        X = F.clauses{i}.data;X = X(:);
        AB = getbase(X);
        K.l = prod(size(X));
        variables = getvariables(X);
        if is_interval(i)
            [lb,ub,cand_rows] = findulb_interval(AB,K);
        else
            [lb,ub,cand_rows] = findulb(AB,K);            
        end
        LU(variables,1) = max([lb LU(variables,1)]')';
        LU(variables,2) = min([ub LU(variables,2)]')';  
        
    elseif F.clauses{i}.type == 3 && nargin==1
        % FIX : Extract from equalities and binary constraints
        X = F.clauses{i}.data;
        AB = getbase(X);
        if size(AB,1)==1 && (AB(1) > 0) && all(AB(2:end)<0)
            % Save this simplex sum a_ix_i == b for later
            % and extract when as many variables as possible have been
            % fixed. Common in portfolio models
            simplex = [simplex i];
        end
        AB = [AB;-AB];
        K.l = 2*prod(size(X));
        variables = getvariables(X);
        if is_interval(i)
            [lb,ub,cand_rows] = findulb_interval(AB,K);
        else
            [lb,ub,cand_rows] = findulb(AB,K);            
        end
        LU(variables,1) = max([lb LU(variables,1)]')';
        LU(variables,2) = min([ub LU(variables,2)]')';
    end
end

semicont = yalmip('semicontvariables');
if ~isempty(semicont)
    for i = 1:length(semicont)
        if LU(semicont(i),1) > 0
            LU(semicont(i),1) = 0;
        end
        if LU(semicont(i),2) < 0
            LU(semicont(i),2) = 0;
        end
    end
end

% Go through the simplex models with as much bounded as possible
for i = simplex
    X = F.clauses{i}.data;
    AB = getbase(X);
    variables = getvariables(X);
    if all(LU(variables,1)>=0)
        simplex_bound = AB(1)./(-AB(2:end));
        LU(variables,2) = min([simplex_bound(:) LU(variables,2)]')';
    end
end

% Try to bound some nonlinear terms
% FIX: complete code
[mt,variable_type] = yalmip('monomtable');
quadratic = find(variable_type == 2);
if ~isempty(quadratic)
    M = mt(quadratic,:);
    [ii,jj,kk] = find(M);
    LU(quadratic,1) = max([LU(quadratic,1) zeros(length(quadratic),1)],[],2);    
    LU(quadratic,2) = max([LU(quadratic,1).^2 LU(quadratic,2).^2],[],2);
end