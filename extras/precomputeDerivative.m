function precompute = precomputeDerivative(model,requested)

Z = double(model.monomtable | model.monomtable);
precompute = cell(length(model.linearindicies),length(model.evaluation_scheme));
%if any(model.deppattern(requested,model.linearindicies(variable)))
M = model.deppattern(requested,:);
active = any(M(:,model.linearindicies),1);
range = 1:length(model.linearindicies);
for variable = range(active(range))% = range
    %dx = dx0;%zeros(length(model.c),1);
    %dx(model.linearindicies(variable)) = 1;
    %if any(model.deppattern(requested,model.linearindicies(variable)))
  %  if active(variable)%any(M(:,model.linearindicies(variable)))
        dx = full(sparse(model.linearindicies(variable),1,1,length(model.c),1));
        for i = 1:length(model.evaluation_scheme)
            precompute{variable,i}=0;
            switch model.evaluation_scheme{i}.group
                case 'eval'
                    for j = model.evaluation_scheme{i}.variables
                        k = model.evalMap{j}.variableIndex;
                        if any(dx(k))
                            if any(requested(model.evalMap{j}.computes))
                                precompute{variable,i}(j)=1;
                              %  if any(requested(model.evalMap{j}.computes))
                                    these = model.evalMap{j}.computes;
                                    dx(these) = 1;
                             %   end
                            end
                        end
                    end
                case 'monom'
                    computed = model.monomials(model.evaluation_scheme{i}.variables);                   
                    hh = any(Z(computed,find(dx)),2); %#ok<FNDSB>
						  j = computed(hh);
						  dx(j(requested(j))) = 1;
                    
                otherwise
            end
            precompute{variable,i} = find(precompute{variable,i});
        end
  %  end
end
