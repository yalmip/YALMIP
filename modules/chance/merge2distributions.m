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
        varC = blkdiag(varA,varB);
    else
        varC = [];
    end
    if ~isempty(facA) && ~isempty(facB)
        facC = blkdiag(facA,facB);
    else
        facC = [];
    end
        
	C.distribution.parameters{1} = 'normal';
    C.distribution.parameters{2} = [A.distribution.parameters{2};B.distribution.parameters{2}];
    C.distribution.parameters{3} = varC;
    C.distribution.parameters{4} = facC;        
else
    for k = 2:length(A.distribution.parameters)
        C.distribution.parameters{k} = [A.distribution.parameters{k};B.distribution.parameters{k}];
    end
end