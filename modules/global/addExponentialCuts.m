function pcut = addExponentialCuts(p)

pcut = p;
if length(p.evalMap) > 0
    pcut = emptyNumericalModel;    
    for i = 1:length(p.evalMap)        
        if isequal(p.evalMap{i}.fcn,'exp')            
            % We have the monomial y = exp(x) in the model
            x = p.evalMap{i}.variableIndex;            
            y = p.evalMap{i}.computes;            
            % Construct the EXPCONE model for 1*exp(x/1) <= y
            % z(2)*exp(z(1)/z(2)) <= z(3)
            F_structemp = spalloc(3,length(p.c)+1,3);            
            F_structemp(1,1+x) = 1;
            F_structemp(2,1) = 1;
            F_structemp(3,1+y) = 1;           
            K.f = 0;
            K.l = 0;
            K.s = 0;
            K.e = 1;
            K.q = 0;
            localModel = createNumericalModel(F_structemp,K);
            pcut = mergeNumericalModels(pcut,localModel);
        end        
    end
    pcut = mergeNumericalModels(p,pcut);
end
