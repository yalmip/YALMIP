function LU = getbounds(F)

K.f = 0;
K.l = 0;
L = [];
U = [];
LU =  yalmip('getbounds',1:yalmip('nvars'));
for i = 1:length(F.clauses)
    if F.clauses{i}.type == 2
        X = F.clauses{i}.data;
        AB = getbase(X);
        K.l = prod(size(X));
        variables = getvariables(X);
        [lb,ub,cand_rows] = findulb(AB,K);
        LU(variables,1) = max([lb LU(variables,1)]')';
        LU(variables,2) = min([ub LU(variables,2)]')';
    elseif F.clauses{i}.type == 3
        % FIX : Extract from equalities and binary constraints
    end
end

binary = yalmip('binvariables');
LU(binary,1) = 0;
LU(binary,2) = 1;

% Try to bound some nonlinear terms
% FIX: complete code
[mt,variable_type] = yalmip('monomtable');
quadratic = find(variable_type == 2);
if ~isempty(quadratic)
    M = mt(quadratic,:);
    for i = 1:size(M,1)
        [ii,jj] = find(M(i,:));
        if length(ii) == 1
            LU(quadratic(i),1) = min([0 LU(jj,1)^2]);
            LU(quadratic(i),2) = max([LU(jj,1)^2 LU(jj,2)^2]);              
        end
    end
end
