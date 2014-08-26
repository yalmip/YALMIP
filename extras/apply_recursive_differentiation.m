function dX = apply_recursive_differentiation(model,x,requested,recursivederivativeprecompute)

global newmodel
persistent dX0
persistent index
persistent mtT
persistent monomTablePattern
persistent ss
persistent lastrequested

% Compute all evaluation-based derivatives df(x)
dxi = [];
dxj = [];
dxs = [];
for i = 1:length(model.evaluation_scheme)
    if isequal(model.evaluation_scheme{i}.group,'eval')
        for j = model.evaluation_scheme{i}.variables
            k = model.evalMap{j}.variableIndex;
            if any(requested(model.evalMap{j}.computes))
                derivative = model.evalMap{j}.properties.derivative(x(k));
                z{i,j} = derivative;
                if i == 1
                    dxj = [dxj model.evalMap{j}.variableIndex];
                    if length(derivative)>1 & length(model.evalMap{j}.computes)==1
                        dxi = [dxi repmat(model.evalMap{j}.computes,1,length(derivative))];
                    else
                        dxi = [dxi model.evalMap{j}.computes];
                    end
                    dxs = [dxs;derivative(:)];
                end
            end
        end
    end
end

% Apply chain-rule. This code is horrible

if newmodel || ~isequal(requested,lastrequested);
    % Some precalc save over iterations
    ss = any(model.deppattern(requested,model.linearindicies),1);
    monomTablePattern = model.monomtable | model.monomtable;
    mtT = model.monomtable';
    [~,dxj] = ismember(dxj,model.linearindicies);
    dxi1 = [dxi model.linearindicies];
    dxj1 = [dxj 1:length(model.linearindicies)];
    dxs1 = [dxs(:)' ones(length(model.linearindicies),1)'];
    dX = sparse(dxi1,dxj1,dxs1,length(model.c),length(model.linearindicies));    
    dX0 = dX;
    index = sub2ind([length(model.c),length(model.linearindicies)],dxi,dxj);   
    newmodel = 0;
    lastrequested = requested;
else
    dX = dX0;
    dX(index) = dxs;
end

newMonoms = [];
for i = 1:length(model.evaluation_scheme)
    switch model.evaluation_scheme{i}.group
        case 'eval'
            if i>1
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
            end
        case 'monom'
            computed = model.monomials(model.evaluation_scheme{i}.variables);
            
            % quadratic and bilinear expressions in the bottom layer of
            % the computational tree can be differentiated very easily
            if i == 1
                BilinearIndex = find(model.variabletype(computed)==1);
                Bilinears = computed(BilinearIndex);
                if ~isempty(Bilinears)
                    x1 = model.BilinearsList(Bilinears,1);
                    x2 = model.BilinearsList(Bilinears,2);
                    dX(sub2ind(size(dX),[Bilinears(:);Bilinears(:)],[x1;x2]))=x([x2;x1]);
                end
                QuadraticIndex = find(model.variabletype(computed)==2);
                Quadratics = computed(QuadraticIndex);
                if ~isempty(Quadratics)
                    x1 =  model.QuadraticsList(Quadratics,1);
                    dX(sub2ind(size(dX),Quadratics',x1))=2*x(x1);
                end
                computed([BilinearIndex QuadraticIndex])=[];
            end
            
            if i > 1
                % We might have inner derivatives from earlier
                isBilinear = model.variabletype==1;
                if all(isBilinear(computed))
                    % Special case quicker. Only bilinear monomials
                    b1 = model.BilinearsList(computed,1);
                    b2 = model.BilinearsList(computed,2);
                    dp = repmat(x(b2),1,length(model.linearindicies)).*dX(b1,:)+repmat(x(b1),1,length(model.linearindicies)).*dX(b2,:);
                    dX(computed,:) = dp;
                else
                    for variable = 1:length(model.linearindicies)
                        if ss(variable)
                            dx = dX(:,variable)';
                            fdX = find(dx);
                            hh = sum(monomTablePattern(computed,fdX),2);
                            for j = computed(find(hh))
                                if requested(j)
                                    if isBilinear(j)
                                        b1 = model.BilinearsList(j,1);
                                        b2 = model.BilinearsList(j,2);
                                        dp = x(b2)*dX(b1,variable)+x(b1)*dX(b2,variable);
                                    else
                                        dp = 0;
                                        monomsj = mtT(:,j)';
                                        active = find(monomsj & dX(:,variable)');
                                        for k = active
                                            monoms = monomsj;
                                            r = monoms(k);
                                            monoms(k) = r-1;
                                            s = find(monoms);
                                            ztemp = prod((x(s)').^(monoms(s)));
                                            aux = dX(k,variable);
                                            dp = dp + r*aux*ztemp;
                                        end
                                    end
                                    dX(j,variable) = real(dp);
                                end
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
