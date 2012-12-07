function Z = jacobiansparsityfromnonlinear(model)

    m = length(model.lb);    
    allA=[model.Anonlinineq; model.Anonlineq];
    jacobianstructure = spalloc(size(allA,1),m,0);    
    depends = allA | allA;   
    for i = 1:size(depends,1)
        vars = find(depends(i,:));
        [ii,vars] = find(model.deppattern(vars,:));
        vars = unique(vars);
        s = size(jacobianstructure,1);
        for j = 1:length(vars)            
            jacobianstructure(i,find(vars(j) == model.linearindicies)) = 1; 
        end      
    end
    allA=[model.A; model.Aeq];
    depends = allA | allA;
    jacobianstructure = [jacobianstructure;depends];
    
    Z = sparse(jacobianstructure);

