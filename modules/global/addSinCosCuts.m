function pcut = addSinCosCuts(p)

pcut = p;
if ~isempty(p.evalMap) 
    sin_ = [];
    cos_ = [];
    for i = 1:length(p.evalMap)
        if isequal(p.evalMap{i}.fcn,'sin')
            sin_ = [sin_;p.evalMap{i}.variableIndex p.evalMap{i}.computes];
        elseif isequal(p.evalMap{i}.fcn,'cos')
            cos_ = [cos_;p.evalMap{i}.variableIndex p.evalMap{i}.computes];
        end
    end
    if ~isempty(sin_) && ~isempty(cos_) && ~isempty(intersect(sin_(:,1),cos_(:,1)))
        for i = 1:size(sin_,1)
            j = find(sin_(i,1)==cos_(:,1));
            if ~isempty(j)
                k1 = sin_(i,2);
                k2 = cos_(j,2);
                k3 = sin_(i,1);
                % x_k1 and x_k2 rpresent sin and cos of x_k3
                                                                 
                % Construct the SOCP model for 1 >= x_k1^2 + x_k2^2
                % i.e. norm([x_k1;x_k2]) <= 1
                % i.e. [1;x_k1;x_k2] in socp cone
                if 0
                    F_structemp = spalloc(3,length(p.c)+1,5);
                    F_structemp(:,1) = [1;0;0];
                    F_structemp(2,1+k1) = 1;
                    F_structemp(3,1+k2) = 1;
                    K.f = 0;
                    K.l = 0;
                    K.s = 0;
                    K.e = 0;
                    K.q = 3;
                    localModel1 = createNumericalModel(F_structemp,K);
                    pcut = mergeNumericalModels(pcut,localModel1);    
                end
                
                % We also add simple linear cuts, in case the lower bound
                % solver doesn't support socps. FIX ME: Decide globally,
                % but make sure it is consistent with stand-alone envelope
                % y <= sqrt(2)-x, y>=-sqrt(2)+x, 
                % y <= sqrt(2)+x, y>=-sqrt(2)-x
                % sqrt(2)-x-y, sqrt(2)+x-y, sqrt(2)-x+y, sqrt(2)+x+y
                F_structemp = spalloc(4,length(p.c)+1,5);
                F_structemp(:,1) = sqrt(2)*[1;1;1;1];
                F_structemp(1,1+k1) = -1;
                F_structemp(1,1+k2) = -1;
                F_structemp(2,1+k1) = 1;
                F_structemp(2,1+k2) = -1;
                F_structemp(3,1+k1) = -1;
                F_structemp(3,1+k2) = 1;
                F_structemp(4,1+k1) = 1;
                F_structemp(4,1+k2) = 1;
                K.f = 0;
                K.l = 4;
                K.s = 0;
                K.e = 0;
                K.q = 0;
                localModel2 = createNumericalModel(F_structemp,K);                               
                pcut = mergeNumericalModels(pcut,localModel2);    
            end
        end
    end    
end