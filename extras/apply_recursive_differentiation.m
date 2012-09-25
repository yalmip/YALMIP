function dX = apply_recursive_differentiation(model,x,requested,recursivederivativeprecompute);

% Compute all evaluation-based derivatives df(x)
for i = 1:length(model.evaluation_scheme)
    if isequal(model.evaluation_scheme{i}.group,'eval')
        for j = model.evaluation_scheme{i}.variables
            k = model.evalMap{j}.variableIndex;
            if any(requested(model.evalMap{j}.computes))
                z{i,j} = model.evalMap{j}.properties.derivative(x(k));
            end
        end
    end
end

% Apply chain-rule. This code is horrible
ss = model.deppattern(requested,model.linearindicies);
ss = ss | ss;
ss = sum(ss,1)>1;

monomTablePattern = sparse(double(model.monomtable | model.monomtable));
mtT = model.monomtable';
if 1 % USE NEW
    dX = sparse(model.linearindicies,1:length(model.linearindicies),ones(length(model.linearindicies),1),length(model.c),length(model.linearindicies)); 
    newMonoms = [];
    for i = 1:length(model.evaluation_scheme)
        switch model.evaluation_scheme{i}.group
            case 'eval'
                for variable = 1:length(model.linearindicies)
                    if ss(variable)
                        for j = recursivederivativeprecompute{variable,i}
                            k = model.evalMap{j}.variableIndex;
                            if length(model.evalMap{j}.computes) == 1
                                dX(model.evalMap{j}.computes,variable) = dX(k,variable)'*z{i,j};
                            else
                                dX(model.evalMap{j}.computes,variable) = dX(k,variable).*z{i,j};
                            end
                        end
                    end
                end
            case 'monom'
                computed = model.monomials(model.evaluation_scheme{i}.variables);
                
                % quadratic and bilinear expressions in the bottom layer of
                % the computational tree can be differentiated very easily 
                if i == 1
                    Bilinears = computed(find(model.variabletype(computed)==1));
                    if ~isempty(Bilinears)
                        x1 = model.BilinearsList(Bilinears,1);
                        x2 = model.BilinearsList(Bilinears,2);
                        dX(sub2ind(size(dX),Bilinears',x1))=x(x2);
                        dX(sub2ind(size(dX),Bilinears',x2))=x(x1);
                    end
                    Quadratics = computed(find(model.variabletype(computed)==2));
                    if ~isempty(Quadratics)
                        x1 =  model.QuadraticsList(Quadratics,1);
                        dX(sub2ind(size(dX),Quadratics',x1))=2*x(x1);
                    end                    
                    Alreadydone = [Bilinears Quadratics];
                else
                    Alreadydone=[];
                end
                
                % These are more complicated than just bilinear, so they
                % have to be taken care of separately
%                 computed1 = setdiff(computed,Alreadydone);   
%                 computed2 = setdiff1D(computed,Alreadydone);   
%                 if ~isequal(computed1,computed2)
%                     1
%                 end
               computed = setdiff(computed,Alreadydone);    
                if i > 1
                    % We might have inner derivatives from earlier
                    for variable = 1:length(model.linearindicies)
                        if ss(variable)
                            dx = dX(:,variable)';
                            fdX = find(dx);
                            hh = sum(monomTablePattern(computed,fdX),2);
                            for j = computed(find(hh))
                                if requested(j)
                                    dp = 0;
                                    monomsj = mtT(:,j)';                                    
                                    for k = find(monomsj & dX(:,variable)')                                                                                
                                        monoms = monomsj;
                                        r = monoms(k);
                                        monoms(k) = r-1;
                                        s = find(monoms);
                                        z = prod((x(s)').^(monoms(s)));                                       
                                        dp = dp + r*dX(k,variable)*z;
                                    end
                                    dX(j,variable) = real(dp);
                                end
                            end
                        end
                    end
                else
                    % We are at the bottom layer, so no inner derivatives
                    for variable = 1:length(model.linearindicies)
                        if ss(variable)
                            dx = dX(:,variable)';
                            fdX = find(dx);
                            hh = sum(monomTablePattern(computed,fdX),2);
                            for j = computed(find(hh))
                                if requested(j)                                  
                                    monomsj = mtT(:,j)';                                    
                                
                                    k = model.linearindicies(variable);
                                    monoms = monomsj;
                                    r = monoms(k);
                                    monoms(k) = r-1;
                                
                                    s = find(monoms);
                                    dp = r*prod((x(s)').^monoms(s));                                    
                                    dX(j,variable) = real(dp);
                                end
                            end
                        end
                    end
                end
                
            otherwise
        end
    end
else
    for variable = 1:length(model.linearindicies)
        dx = zeros(length(model.c),1);
        dx(model.linearindicies(variable)) = 1;
        if ss(variable)
            for i = 1:length(model.evaluation_scheme)
                switch model.evaluation_scheme{i}.group
                    case 'eval'
                        for j = recursivederivativeprecompute{variable,i}
                            k = model.evalMap{j}.variableIndex;
                            if length(model.evalMap{j}.computes) == 1
                                dx(model.evalMap{j}.computes) = dx(k)'*z{i,j};
                            else
                                dx(model.evalMap{j}.computes) = dx(k).*z{i,j};
                            end
                        end
                    case 'monom'
                        computed = model.monomials(model.evaluation_scheme{i}.variables);
                        % hh = model.monomtable(computed,:);
                        % hh = double(hh | hh)*(dx | dx);
                        hh = sum(monomTablePattern(computed,find(dx)),2);
                        for j = computed(find(hh))
                            if requested(j)
                                dp = 0;
                                monomsj = model.monomtable(j,:);
                                for k = find(dx' & monomsj)
                                    monoms = monomsj;
                                    monoms(k) = 0;
                                    r = model.monomtable(j,k);
                                    s = find(monoms);
                                    dp = dp + r*x(k)^(r-1)*dx(k)*prod((x(s)').^monoms(s));
                                end
                                dx(j) = real(dp);
                            end
                        end
                        
                    otherwise
                end
            end
        end
        %dX = [dX dx];
        dX(:,variable) = dx;
    end
end