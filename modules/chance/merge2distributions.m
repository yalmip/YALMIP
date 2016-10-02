function C = merge2distributions(A,B)

C = A;
C.variables = [A.variables;B.variables];
if any(strcmp(A.distribution.parameters{1},{'normal','normalf','normalm'})) || any(strcmp(A.distribution.parameters{1},{'normal','normalf','normalm'}))
    
    % A bit messy as there are two additional forms for normal
    % distributions, used to define mv normals, and normals with a
    % factorized covariance
    
    % First, normalize to matrix format
    dimAvar = size(A.distribution.parameters{3});
    dimBvar = size(B.distribution.parameters{3});
    if dimAvar(1) ~= dimAvar(2)
        varA = diag(A.distribution.parameters{3});
    else
        varA = A.distribution.parameters{3};
        if numel(varA)==1 && numel(A.variables)>1
            varA = diag(repmat(varA,numel(A.variables),1));
        end
    end
    if dimBvar(1) ~= dimBvar(2)
        varB = diag(B.distribution.parameters{3});
    else
        varB = B.distribution.parameters{3};
        if numel(varB)==1 && numel(B.variables)
            varB = diag(repmat(varB,numel(A.variables),1));
        end
    end
    % Is A in factor form, but not B, and vice versa. If so, put the other
    % one in factor form too
    if strcmp(A.distribution.parameters{1},'normalf') &&  any(strcmp(B.distribution.parameters{1},{'normal','normalm'}))
        varB = chol(varB);
    elseif strcmp(B.distribution.parameters{1},'normalf') &&  any(strcmp(A.distribution.parameters{1},{'normal','normalm'}))
        varA = chol(varA);
    end
    varC = blkdiag(varA,varB);
    
    % Sort out various combinations of normal stuff in different forms
    if strcmp(A.distribution.parameters{1},'normal') && strcmp(B.distribution.parameters{1},'normal')
        varC = diag(varC);
        C.distribution.parameters{1} = 'normal';
        C.distribution.parameters{2} = [A.distribution.parameters{2};B.distribution.parameters{2}];
        C.distribution.parameters{3} = diag(varC);
    elseif strcmp(A.distribution.parameters{1},'normalf') ||  any(strcmp(B.distribution.parameters{1},'normalf'))
        C.distribution.parameters{1} = 'normalf';
        C.distribution.parameters{2} = [A.distribution.parameters{2};B.distribution.parameters{2}];
        C.distribution.parameters{3} = varC;
    else
        C.distribution.parameters{1} = 'normalm';
        C.distribution.parameters{2} = [A.distribution.parameters{2};B.distribution.parameters{2}];
        C.distribution.parameters{3} = varC;
    end
else
    for k = 2:length(A.distribution.parameters)
        C.distribution.parameters{k} = [A.distribution.parameters{k};B.distribution.parameters{k}];
    end
end