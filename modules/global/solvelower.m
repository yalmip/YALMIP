function [output,cost,psave,timing] = solvelower(p,options,lowersolver,xmin,upper,timing)

psave = p;
removeThese = find(p.InequalityConstraintState==inf);
p.F_struc(p.K.f + removeThese,:) = [];
p.K.l = p.K.l - length(removeThese);

removeThese = find(p.EqualityConstraintState==inf);
p.F_struc(removeThese,:) = [];
p.K.f = p.K.f - length(removeThese);

p_cut = p;

% Add cuts for x^2==x and variants for box QP
p_cut.F_struc = [p_cut.concavityEqualities;p_cut.F_struc];
p_cut.K.f = p_cut.K.f + size(p_cut.concavityEqualities,1);

if p.options.bmibnb.cut.bilinear
    p_cut = addBilinearVariableCuts(p_cut);
end

if p.options.bmibnb.cut.normbound
  	p_cut = addNormBoundCut(p_cut);
end
if p.options.bmibnb.cut.evalvariable
    p_cut = addEvalVariableCuts(p_cut);
    psave.evalMap = p_cut.evalMap;
end
if p.options.bmibnb.cut.monomial
	p_cut = addMonomialCuts(p_cut);
end
if p.options.bmibnb.cut.complementarity
  	p_cut = addComplementarityCuts(p_cut);
end
if p.solver.lowersolver.constraint.inequalities.secondordercone.linear && p.options.bmibnb.cut.quadratic
    p_cut = addQuadraticCuts(p_cut);
end
if p.options.bmibnb.cut.exponential
    p_cut = addExponentialCuts(p_cut);
end
if p.options.bmibnb.cut.sincos
    p_cut = addSinCosCuts(p_cut);
end
if p.options.bmibnb.cut.monomialtower
    p_cut = addMonomialTowerCuts(p_cut);
end
if ~isempty(p_cut.binary_variables) || ~isempty(p_cut.integer_variables)
    if ~isempty(p_cut.K.s) & p_cut.K.s(1) > 0
        if isequal(p_cut.solver.lowercall,'callmosek')
            % Mosek SDP module does not support binary
            p_cut.binary_variables = [];
            p_cut.integer_variables = [];
            p_cut.binary_variables = [];
        end
    end
end
% **************************************
% SOLVE NODE PROBLEM
% **************************************
if any(p_cut.ub+1e-8<p_cut.lb)
    output.problem=1;
    cost = inf;
else
    % We are solving relaxed problem (penbmi might be local solver)
    p_cut.monomtable = eye(length(p_cut.c));
    
    % If lower solver can handle convex quadratic, we should exploit
    if p.solver.lowersolver.objective.quadratic.convex && ~isempty(p.shiftedQP)
        % Use pre-computed SDP-shifted model
        p_cut.Q = p.shiftedQP.Q;
        p_cut.c = p.shiftedQP.c;       
        p_cut.f = p.shiftedQP.f;          
    elseif p.solver.lowersolver.objective.quadratic.convex 
        % If we have a QP solver, we can at least try to strengthen the
        % relaxation by using the positive diagonal terms
        [p_cut.Q, p_cut.c] =  compileQuadratic(p_cut.c,p,3);
    end
    
    fixed = p_cut.lb >= p_cut.ub;
    if nnz(fixed) == length(p.c)
        % All variables are fixed to a bound
        output.Primal = p.lb;
        res = constraint_residuals(p,output.Primal);
        eq_ok = all(res(1:p.K.f)>=-p.options.bmibnb.eqtol);
        iq_ok = all(res(1+p.K.f:end)>=-p.options.bmibnb.pdtol);
        feasible = eq_ok & iq_ok;
        if feasible
            output.problem = 0;
        else
            output.problem = 1;
        end
        cost = output.Primal'*p.Q*output.Primal + p.c'*output.Primal + p.f;
    else
        
        if nnz(fixed)==0
            
            if ~isempty(p_cut.bilinears) & 0
                top = size(p_cut.F_struc,1);
                if length(p_cut.K.s)==1 & p_cut.K.s(1)==0
                    p_cut.K.s = [];
                end
                usedterms = zeros(size(p_cut.bilinears,1),1);
                for i = 1:size(p_cut.bilinears,1)
                    if ~usedterms(i)
                        windex = p_cut.bilinears(i,1);
                        xindex = p_cut.bilinears(i,2);
                        yindex = p_cut.bilinears(i,3);
                        if xindex ~=yindex
                            % OK, we have a bilinear term
                            xsquaredindex = find(p_cut.bilinears(:,2)==xindex & p_cut.bilinears(:,3)==xindex);
                            ysquaredindex = find(p_cut.bilinears(:,2)==yindex & p_cut.bilinears(:,3)==yindex);
                            if ~isempty(xsquaredindex) & ~isempty(ysquaredindex)
                                usedterms(i) = 1;
                                usedterms(xsquaredindex) = 1;
                                usedterms(ysquaredindex) = 1;
                                xsquaredindex =  p_cut.bilinears(xsquaredindex,1);
                                ysquaredindex =  p_cut.bilinears(ysquaredindex,1);
                                if 0
                                    Z = zeros(9,size(p_cut.F_struc,2));
                                    Z(1,xsquaredindex+1) = 1;
                                    Z(2,windex+1) = 1;
                                    Z(4,windex+1) = 1;
                                    Z(5,ysquaredindex+1) = 1;
                                    Z(3,xindex+1) = 1;
                                    Z(7,xindex+1) = 1;
                                    Z(6,yindex+1) = 1;
                                    Z(8,yindex+1) = 1;
                                    Z(9,1)=1;
                                else
                                    xL = p.lb(xindex);
                                    yL = p.lb(yindex);
                                    
                                    Z = zeros(9,size(p_cut.F_struc,2));
                                    Z(1,xsquaredindex+1) = 1;
                                    Z(2,windex+1) = 1;
                                    Z(4,windex+1) = 1;
                                    Z(5,ysquaredindex+1) = 1;
                                    Z(3,xindex+1) = 1;
                                    Z(7,xindex+1) = 1;
                                    Z(6,yindex+1) = 1;
                                    Z(8,yindex+1) = 1;
                                    Z(9,1)=1;
                                    Z(3,1) = -xL;
                                    Z(7,1) = -xL;
                                    Z(6,1) = -yL;
                                    Z(8,1) = -yL;
                                    
                                    Z(1,xindex+1) = -2*xL;
                                    Z(5,yindex+1) = -2*yL;
                                    
                                    Z(1,1) = xL^2;
                                    Z(5,1) = yL^2;
                                    
                                    Z(4,xindex+1) = -yL;
                                    Z(4,yindex+1) = -xL;
                                    Z(4,1) = xL*yL;
                                    
                                    Z(2,xindex+1) = -yL;
                                    Z(2,yindex+1) = -xL;
                                    Z(2,1) = xL*yL;
                                    
                                    
                                end
                                p_cut.F_struc = [p_cut.F_struc;Z];
                                p_cut.K.s = [p_cut.K.s 3];
                            end
                        end
                    end
                end
            end
            
            p_cut.linearindicies = 1:length(p.c);
            p_cut.nonlinearindicies = [];
            p_cut.variabletype = zeros(1,length(p.c));
            p_cut.monomtable = speye(length(p.c));
            p_cut.deppattern = eye(length(p.c));
            p_cut.linears = 1:length(p.c);
            p_cut.bilinears = [];
            p_cut.nonlinears = [];
            p_cut.monomials = [];
            p_cut.evaluation_scheme = [];
            
            tstart = tic;     
            p_cut = pruneUnsupportedCuts(p_cut);            
            output = feval(lowersolver,removenonlinearity(p_cut));
            psave.counter.lowersolved = psave.counter.lowersolved + 1;
            timing.lowersolve = timing.lowersolve + toc(tstart);
            if length(output.Primal) == length(p_cut.c)
                cost = output.Primal'*p_cut.Q*output.Primal + p_cut.c'*output.Primal + p_cut.f;
                % Minor clean-up
                pp=p;
                output.Primal(output.Primal<p.lb) = p.lb(output.Primal<p.lb);
                output.Primal(output.Primal>p.ub) = p.ub(output.Primal>p.ub);
                x=output.Primal;
            else
                cost = nan;
                x = nan(length(p_cut.c),1);
                output.Primal = x;
            end
            return
        else
            pp = p_cut;
            removethese = fixed;
            if ~isempty(p_cut.F_struc)
                p_cut.F_struc(:,1)=p_cut.F_struc(:,1)+p_cut.F_struc(:,1+find(fixed))*p_cut.lb(fixed);
                p_cut.F_struc(:,1+find(fixed))=[];
                
                rf = find(~any(p_cut.F_struc,2));
                rf = rf(rf<=(p_cut.K.f + p_cut.K.l));
                p_cut.F_struc(rf,:) = [];
                p_cut.K.l = p_cut.K.l - nnz(rf>p_cut.K.f);
                p_cut.K.f = p_cut.K.f - nnz(rf<=p_cut.K.f);
            elseif size(p_cut.F_struc,2)>0
                p_cut.F_struc(:,1+find(fixed))=[]; % can be 0x(n+1)
            end
            p_cut.c(removethese)=[];
            if nnz(p_cut.Q)>0
                p_cut.c = p_cut.c + 2*p_cut.Q(find(~removethese),find(removethese))*p_cut.lb(removethese);
                p_cut.Q(:,find(removethese))=[];
                p_cut.Q(find(removethese),:)=[];
            else
                p_cut.Q = spalloc(length(p_cut.c),length(p_cut.c),0);
            end
            
            if ~isempty(p_cut.binary_variables)
                new_bin = [];
                new_var = find(~fixed);
                for i = 1:length(p_cut.binary_variables)
                    temp = find(p_cut.binary_variables(i) == new_var);
                    new_bin =  [new_bin temp(:)'];
                end
                p_cut.binary_variables = new_bin;
            end
            if ~isempty(p_cut.integer_variables)
                new_bin = [];
                new_var = find(~fixed);
                for i = 1:length(p_cut.integer_variables)
                    temp = find(p_cut.integer_variables(i) == new_var);
                    new_bin =  [new_bin temp(:)'];
                end
                p_cut.integer_variables = new_bin;
            end
            
            p_cut.lb(removethese)=[];
            p_cut.ub(removethese)=[];
            p_cut.x0(removethese)=[];
            p_cut.monomtable(:,find(removethese))=[];
            p_cut.monomtable(find(removethese),:)=[];
            p_cut.variabletype(removethese) = [];
            
            if p_cut.solver.lowersolver.constraint.integer == 0 && ~(isempty(p_cut.binary_variables) && isempty(p_cut.integer_variables))
                p_cut.integer_variables = [];
                p_cut.binary_variables = [];
            end
            
            % The model can become absolutely trivial in some case
            % For instance in ex9_2_2 everything is presolved
            if nnz(p_cut.c)==0 & nnz(p_cut.Q)==0 & size(p_cut.F_struc,1)==0
                % No objective and no constraints
                if all(p_cut.lb <= p_cut.ub)
                    output.Primal = (p_cut.lb + p_cut.ub)/2;
                    output.Primal(isinf(p_cut.lb) & isinf(p_cut.ub))=0;
                    output.problem = 0;
                else
                    output.Primal = zeros(length(p_cut.lb),1);
                    output.problem = 1;
                end
            else
                try
                    tstart = tic;                    
                    output = feval(lowersolver,removenonlinearity(p_cut));
                    psave.counter.lowersolved = psave.counter.lowersolved + 1;
                    timing.lowersolve = timing.lowersolve + toc(tstart);
                catch
                    1
                end
            end
            x=full(p.c*0);
            x(removethese)=p.lb(removethese);
            x(~removethese)=output.Primal;
            output.Primal = x;
            cost = output.Primal'*pp.Q*output.Primal + pp.c'*output.Primal + p.f;
        end
    end
end

function p = pruneUnsupportedCuts(p)

if p.K.s > 0 & ~p.solver.lowersolver.constraint.inequalities.semidefinite.linear
    % Remove SDP cuts
    error('FIXME')
end

if p.K.e > 0 & ~p.solver.lowersolver.exponentialcone
    % Remove EXPCONE cuts
    p.F_struc(end-p.K.e*3+1:1:end,:) = [];    
    p.K.e = 0;
end
