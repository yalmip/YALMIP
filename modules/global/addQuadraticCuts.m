function pcut = addQuadraticCuts(p)

pcut = p;
if length(p.Quadratics) > 0
    pcut = emptyNumericalModel;
    
    for y = p.Quadratics(1:length(p.Quadratics))
        
        % We have the monomial y = x^2 in the model
        x = p.QuadraticsList(y,1);
        
        % Construct the SOCP model for y >= x^2
        % i.e. norm([2*x;1-y]) <= 1 + y
        % i.e. [1 + y;2x;1-y] in socp cone
        
        F_structemp = spalloc(3,length(p.c)+1,5);
        F_structemp(:,1) = [1;0;1];
        F_structemp(2,1+x) = 2;
        F_structemp(1,1+y) = 1;
        F_structemp(3,1+y) = -1;
        K.f = 0;
        K.l = 0;        
        K.s = 0;
        K.e = 0;
        K.q = 3;
        localModel = createNumericalModel(F_structemp,K);
        pcut = mergeNumericalModels(pcut,localModel);                   
    end    
    pcut = mergeNumericalModels(p,pcut);
end

