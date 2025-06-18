function C = merge2distributions(A,B)

C = A;
C.variables = [A.variables;B.variables];
normalvariants = {'normal','mvnrnd','mvrndfactor'};
if any(strcmp(A.distribution.parameters{1},normalvariants)) || any(strcmp(A.distribution.parameters{1},normalvariants))
    
    % A bit messy as univariate and multivariate are defined differently in
    % statistics toolbox and then we have a YALMIP specific also 
    % normal, scalar elementwise normal (parameterized in std. dev!)  
    % mvrnd, multivariable normal, (parameterized in covariance!)
    % mvrndfactor, m.v. normal, (parameterized in factor covariance)
    
    % To deal with this, we will have three paramters related to Gaussians
    % par{1} mean
    % par{2} covariance
    % par{3} standard deviation/factored covariance
          
    varA = A.distribution.parameters{3};
    varB = B.distribution.parameters{3};
    facA = A.distribution.parameters{4};
    facB = B.distribution.parameters{4};        
    if ~isempty(varA) && ~isempty(varB)
        if isa(varA,'cell')
            for j = 1:length(varA)
                varC{j} = blkdiag(varA{j},varB{j});
            end
        else
            varC = blkdiag(varA,varB);
        end
    else
        varC = [];
    end
    if ~isempty(facA) && ~isempty(facB)
        if isa(facA,'cell')
            for j = 1:length(varA)
                facC{j} = blkdiag(facA{j},facB{j});
            end
        else
            facC = blkdiag(facA,facB);
        end
    else
        facC = [];
    end
    meanA = A.distribution.parameters{2};
    meanB = B.distribution.parameters{2};
    if isa(meanA,'cell')
        for j = 1:length(meanA)
                meanC{j} = [meanA{j};meanB{j}];;
        end 
    else
        meanC = [meanA;meanB];
    end
        
	C.distribution.parameters{1} = 'normal';
    C.distribution.parameters{2} = meanC;
    C.distribution.parameters{3} = varC;
    C.distribution.parameters{4} = facC;        
else
    for k = 2:length(A.distribution.parameters)
        C.distribution.parameters{k} = [A.distribution.parameters{k};B.distribution.parameters{k}];
    end
end