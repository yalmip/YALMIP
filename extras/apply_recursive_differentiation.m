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
    if strcmp(model.evaluation_scheme{i}.group,'eval')
        for j = model.evaluation_scheme{i}.variables
            if any(requested(model.evalMap{j}.computes))
                k = model.evalMap{j}.variableIndex;
                derivative = model.evalMap{j}.properties.derivative(x(k));
                z{i,j} = derivative;
                if i == 1
                    dxj = [dxj model.evalMap{j}.variableIndex];
                    if length(derivative)>1 && length(model.evalMap{j}.computes)==1
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
    monomTablePattern = logical(model.monomtable);
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

for i = 1:length(model.evaluation_scheme)
	group = model.evaluation_scheme{i}.group;
	if strcmp(group, 'eval')
		if i>1
			for variable = find(ss)
				for j = recursivederivativeprecompute{variable,i}
					 k = model.evalMap{j}.variableIndex;
					 if length(model.evalMap{j}.computes) == 1
						  dX(model.evalMap{j}.computes,variable) = dX(k,variable)'*z{i,j};
					 else
					  val = dX(k,variable).*z{i,j};
					  [idx,~,vdx] = find(val);
					  if ~isempty(vdx)
							dX(model.evalMap{j}.computes(idx),variable) = vdx;
					  end
					 end
				end
			end
		end
	elseif strcmp(group, 'monom')
		computed = model.monomials(model.evaluation_scheme{i}.variables);

		% quadratic and bilinear expressions in the bottom layer of
		% the computational tree can be differentiated very easily
		if i == 1
			 BilinearIndex = find(model.variabletype(computed)==1);
			 Bilinears = computed(BilinearIndex);
			 if ~isempty(Bilinears)
				  x1 = model.BilinearsList(Bilinears,1);
				  x2 = model.BilinearsList(Bilinears,2);
				  [~,x1loc] = ismember(x1,model.linearindicies);
				  [~,x2loc] = ismember(x2,model.linearindicies);
				  dX(sub2ind(size(dX),[Bilinears(:);Bilinears(:)],[x1loc;x2loc]))=x([x2;x1]);
			 end
			 QuadraticIndex = find(model.variabletype(computed)==2);
			 Quadratics = computed(QuadraticIndex);
			 if ~isempty(Quadratics)
				  x1 =  model.QuadraticsList(Quadratics,1);
				  [~,x1loc] = ismember(x1,model.linearindicies);
				  dX(sub2ind(size(dX),Quadratics',x1loc))=2*x(x1);
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
				  compMTP = monomTablePattern(computed,:);
				  for variable = find(ss)
						[dXRows,~,dXVals] = find(dX(:,variable));
						hh = any(compMTP(:,dXRows),2);
						compInd = computed(hh);
						newIndices = compInd(requested(compInd)); % all values for loop
						% Create a working copy of dX(:,variable). Real size isn't
						% length(dXRows)+length(newIndices) but length(unique([dXRows;newIndices']).
						% But it's more efficient to allocate some more memory than calculate the true length.
						dx = sparse(dXRows,1,dXVals, size(dX,1),1, length(dXRows)+length(newIndices));
						for j = newIndices(:)'
							if isBilinear(j)
								 b1 = model.BilinearsList(j,1);
								 b2 = model.BilinearsList(j,2);
								 dp = x(b2)*full(dx(b1))+x(b1)*full(dx(b2));
							else
								 dp = 0;
								 monomsj = mtT(:,j);
								 active = find(monomsj & dx);
								 for k = active(:)'
									monoms = monomsj;
									[s,~,v] = find(monoms);
									k2s = s==k;
									r = v(k2s);
									v(k2s) = r-1;
									ztemp = prod(x(s).^v);
									aux = full(dx(k));
									dp = dp + r*aux*ztemp;
								 end
							end
							dx(j) = real(dp);
						end
						dX(:,variable) = dx(:);
				  end
			 end
		else
			% We are at the bottom layer, so no inner derivatives
			compMTP = monomTablePattern(computed,:);
			for variable = find(ss)
				dx = dX(:,variable);
				fdX = find(dx); % logical(dx); for sparse-matrices, FIND is about as fast.
				hh = any(compMTP(:,fdX),2); %#ok<FNDSB>
				compInd = computed(hh);
				for j = compInd(requested(compInd))
					k = model.linearindicies(variable);
					monoms = mtT(:,j);
					
					[s,~,v] = find(monoms);
					k2s = s==k;
					r = v(k2s);
					v(k2s) = r-1;
					dp = r*prod(x(s).^v);
					dX(j,variable) = real(dp);
				end
			end
		end
	end
end
end
