function [model,keptvariablesIndex] = eliminatevariables(model,removedparameters,value,parameters)
    
if isempty(removedparameters) 
    keptvariablesIndex = 1:length(model.c);
    return
end

keptvariablesIndex = 1:length(model.c);

rmvmonoms = model.rmvmonoms;
newmonomtable = model.newmonomtable;

% Assign evaluation-based values
if length(model.evalParameters > 0)
    alls = [removedparameters(:)' model.evalParameters(:)'];  
    removedP = [];
    removedE = [];
    for j = 1:length(removedparameters)
        for i = 1:length(model.evalMap)
            if isequal(model.evalMap{i}.variableIndex,removedparameters(j))
                p = value(j);
                index = find(model.evalMap{i}.computes == alls);
                value(index) = feval(model.evalMap{i}.fcn,p);
                removedP = [removedP model.evalMap{i}.computes];
                removedE = [removedE i];            
            end
        end
    end
    model.evalVariables(removedE)=[];
    model.evalMap = {model.evalMap{setdiff(1:length(model.evalMap),removedE)}};
    model.evalParameters = setdiff(model.evalParameters,removedP);
end

% Evaluate fixed monomial terms
aux = model.precalc.aux;
if ~isempty(model.precalc.jj1)
    z = value(model.precalc.jj1).^rmvmonoms(model.precalc.index1);
    aux(model.precalc.index1) = z;
end
% Assign simple linear terms
if ~isempty(model.precalc.jj2)
    aux(model.precalc.index2) = value(model.precalc.jj2);
end

monomvalue = prod(aux,2);

removethese = model.removethese;
keepingthese = model.keepingthese;

value = monomvalue(removethese);
monomgain = monomvalue(keepingthese);

if ~isempty(model.F_struc)
    model.F_struc(:,1) = model.F_struc(:,1)+model.F_struc(:,1+removethese)*value;
    model.F_struc(:,1+removethese) = [];
    model.F_struc = model.F_struc*diag(sparse([1;monomgain]));
end
infeasible = 0;
if model.K.f > 0
    candidates = find(~any(model.F_struc(1:model.K.f,2:end),2));
    %candidates = find(sum(abs(model.F_struc(1:model.K.f,2:end)),2) == 0);
    if ~isempty(candidates)
        % infeasibles = find(model.F_struc(candidates,1)~=0);
        if find(model.F_struc(candidates,1)~=0,1)%;~isempty(infeasibles)
            model.infeasible = 1;
            return
        else
            model.F_struc(candidates,:) = [];
            model.K.f = model.K.f - length(candidates);
        end
    end
end
if model.K.l > 0
    candidates = find(~any(model.F_struc(model.K.f + (1:model.K.l),2:end),2));
    %candidates = find(sum(abs(model.F_struc(model.K.f + (1:model.K.l),2:end)),2) == 0);
    if ~isempty(candidates)
        
        %if find(model.F_struc(model.K.f + candidates,1)<0,1)
        if any(model.F_struc(model.K.f + candidates,1)<-1e-14)
            model.infeasible  = 1;
            return
        else
            model.F_struc(model.K.f + candidates,:) = [];
            model.K.l = model.K.l - length(candidates);
        end
    end
end
if model.K.q(1) > 0
    removeqs = [];
    removeRows = [];
    top = model.K.f + model.K.l + 1;  
    F_struc = model.F_struc(top:top+sum(model.K.q)-1,:);   
    top = 1;
    if all(any(F_struc(:,2:end),2))
        % There is still something on every row in all SOCPs, we don't have
        % to search for silly SOCPS ||0|| <= constant
    else
        for i = 1:length(model.K.q)
            rows = top:top+model.K.q(i)-1;
            v = F_struc(rows,:);
            if nnz(v(:,2:end))==0
                if norm(v(2:end,1)) > v(1,1)
                    model.infeasible = 1;
                    return
                else
                    removeqs = [removeqs;i];
                    removeRows = [removeRows;model.K.f+model.K.l+rows];
                end
            end
            top = top + model.K.q(i);
        end
        model.K.q(removeqs)=[];
        model.F_struc(removeRows,:)=[];
        if isempty(model.K.q)
            model.K.q = 0;
        end
    end
end
if model.K.s(1) > 0  
    % Nonlinear semidefinite program with parameter
    top = model.K.f + model.K.l + sum(model.K.q) + 1;
    removeqs = [];
    removeRows = [];
    for i = 1:length(model.K.s)
        n = model.K.s(i);
        rows = top:top+n^2-1;
        v = model.F_struc(rows,:);
        if nnz(v)==0
            removeqs = [removeqs;i];
            removeRows = [removeRows;rows];
        elseif nnz(v(:,2:end))==0
            Q = reshape(v(:,1),n,n);
            used = find(any(Q));Qred=Q(:,used);Qred = Qred(used,:);
            [~,p] = chol(Qred);
            if p
                model.infeasible = 1;
                return
            else
                removeqs = [removeqs;i];
                removeRows = [removeRows;rows];
            end
        end
        top = top + n^2;
    end
    model.K.s(removeqs)=[];
    model.F_struc(removeRows,:)=[];
    if isempty(model.K.s)
        model.K.s = 0;
    end
end

model.f = model.f + model.c(removethese)'*value;

model.c(removethese)=[];
if nnz(model.Q)>0
    model.c = model.c + 2*model.Q(keepingthese,removethese)*value;
end

if nnz(model.Q)==0
    model.Q = spalloc(length(model.c),length(model.c),0);
else
    model.Q(removethese,:) = [];
    model.Q(:,removethese) = [];
end
model.c = model.c.*monomgain;
keptvariablesIndex(removethese) = [];

model.lb(removethese)=[];
model.ub(removethese)=[];

newmonomtable(:,removethese) = [];
newmonomtable(removethese,:) = [];

if ~isequal(newmonomtable,model.precalc.newmonomtable)%~isempty(removethese)
    skipped = [];
    alreadyAdded = zeros(1,size(newmonomtable,1));      
    %[ii,jj,kk] = unique(newmonomtable*gen_rand_hash(0,size(newmonomtable,2),1),'rows','stable');
    [ii,jj,kk,skipped] = stableunique(newmonomtable*gen_rand_hash(0,size(newmonomtable,2),1));   
    S = sparse(kk,1:length(kk),1);
   % skipped = setdiff(1:length(kk),jj);
    model.precalc.S = S;
    model.precalc.skipped = skipped;
    model.precalc.newmonomtable = newmonomtable;
    model.precalc.blkOneS = blkdiag(1,S');
else
    S = model.precalc.S;
    skipped = model.precalc.skipped;
end
model.c = S*model.c;
%model.F_struc2 = [model.F_struc(:,1) (S*model.F_struc(:,2:end)')'];
if ~isempty(model.F_struc)
    model.F_struc = model.F_struc*model.precalc.blkOneS;%blkdiag(1,S');
end
%norm(model.F_struc-model.F_struc2)

model.lb(skipped) = [];
model.ub(skipped) = [];
newmonomtable(skipped,:) = [];
newmonomtable(:,skipped) = [];


if nnz(model.Q)==0
    model.Q = spalloc(length(model.c),length(model.c),0);
else
    model.Q(:,skipped)=[];
    model.Q(skipped,:)=[];
end

keptvariablesIndex(skipped) = [];

model.monomtable = newmonomtable;
model = compressModel(model);

x0wasempty = isempty(model.x0);
model.x0 = zeros(length(model.c),1);

% Try to reduce to QP
[model,keptvariablesIndex,newmonomtable] = setupQuadratics(model,keptvariablesIndex,newmonomtable);

if nnz(model.Q) > 0
    if  model.solver.objective.quadratic.convex == 0 &  model.solver.objective.quadratic.nonconvex == 0
        error('The objective instantiates as a quadratic after fixing parameters, but this is not directly supported by the solver. YALMIP will not reformulate models if they structurally change in call to optimizer. A typical trick to circumvent this is to define a new set of variable e, use the quadratic function e''e, and add an equality constraint e = something. The SOCP formulation can then be done a-priori by YALMIP.');
    end
end

if x0wasempty
    model.x0 = [];
end

% Remap indicies
if ~isempty(model.integer_variables)
    temp=ismember(keptvariablesIndex,model.integer_variables);
    model.integer_variables = find(temp);  
end
if ~isempty(model.binary_variables)
    temp=ismember(keptvariablesIndex,model.binary_variables);
    model.binary_variables = find(temp);    
end
if ~isempty(model.semicont_variables)
   temp=ismember(keptvariablesIndex,model.semicont_variables);
   model.semicont_variables = find(temp);  
end
if ~isempty(model.aux_variables)
   temp=ismember(keptvariablesIndex,model.aux_variables);
   model.aux_variables = find(temp);  
end

model.used_variables = model.used_variables(keptvariablesIndex);

% Check if there are remaining strange terms. This occurs in #152
% FIXME: Use code above recursively instead...
if ~model.solver.objective.sigmonial & any(model.variabletype == 4)
    % Bugger. There are signomial terms left, despite elimination, and the
    % solver does not handle this. YALMIP has introduced an intermediate
    % variable which is a nasty function of the parameter.      
    signomials = find(model.variabletype == 4);
    involved = [];
    for i = 1:length(signomials)
        m = model.monomtable(signomials(i),:);
        involved = [involved;find(m ~= fix(m) | m < 0)];
    end
    involved = unique(involved);
    [lb,ub] = findulb(model.F_struc,model.K,model.lb,model.ub);
    if all(lb(involved) == ub(involved))
        % Now add equality constraints to enforce       
        for i = signomials
            m = model.monomtable(i,:);
            involved = find(m ~= fix(m) | m < 0);
            gain = lb(involved).^m(involved);
            s = zeros(1,size(model.F_struc,2));
            multiplies = setdiff(find(m),involved);
            s(i+1) = 1;
            s(multiplies+1) = -gain;
            model.F_struc = [s;model.F_struc];
            model.K.f = model.K.f + 1;
        end
    else
        error('Did not manage to instatiate model. Complicating terms remaining');
    end
end

function model = compressModel(model)
model.variabletype = zeros(size(model.monomtable,1),1)';
nonlinear = sum(model.monomtable,2)~=1 | sum(model.monomtable~=0,2)~=1;
if ~isempty(nonlinear)
    model.variabletype(nonlinear) = 3;
    quadratic = sum(model.monomtable,2)==2;
    model.variabletype(quadratic) = 2;
    bilinear = max(model.monomtable,[],2)<=1;
    model.variabletype(bilinear & quadratic) = 1;
    sigmonial = any(0>model.monomtable,2) | any(model.monomtable-fix(model.monomtable),2);
    model.variabletype(sigmonial) = 4;
end




function [model,keptvariables,newmonomtable] = setupQuadratics(model,keptvariables,newmonomtable);
if any(model.variabletype) & all(model.variabletype <= 2)
    monomials = find(model.variabletype);
    if nnz(model.F_struc(:,1+monomials))==0
        if all(isinf(model.lb(monomials)))
            if all(isinf(model.ub(monomials)))     
                % OK, the quadratic/bilinear terms only enter in objective,
                % so let us try to construct Q
                if isempty(model.precalc.Qmap)
                ii = [];
                jj = [];
                kk = [];
                for k = monomials(:)'
                    i = find(model.monomtable(k,:));
                    if model.variabletype(k)==1
                        model.Q(i(1),i(2)) = model.Q(i(1),i(2)) + model.c(k)/2;
                        model.Q(i(2),i(1)) = model.Q(i(1),i(2));
                        ii = [ii;i(1)];
                        jj = [jj;i(2)];                        
                        kk = [kk;k];
                    else
                        model.Q(i,i) = model.Q(i,i) + model.c(k);
                        ii = [ii;i];
                        jj = [jj;i];
                        kk = [kk;k];
                    end
                end
                model.precalc.Qmap.M = sparse(sub2ind(size(model.Q),ii,jj),kk,1,prod(size(model.Q)),length(model.c));
                model.precalc.Qmap.monomials = monomials;
                model.precalc.Qmap.monomtable = model.monomtable;
                else
                    m =  model.precalc.Qmap.M*sparse(model.c);
                    m = reshape(m,sqrt(length(m)),[]);
                    model.Q = (m+m')/2;
                end
                model.c(monomials)=[];
                model.F_struc(:,1+monomials) = [];
                model.lb(monomials) = [];
                model.ub(monomials) = [];
                newmonomtable(monomials,:) = [];
                newmonomtable(:,monomials) = [];
                model.monomtable = newmonomtable;
                model.Q(:,monomials) = [];
                model.Q(monomials,:) = [];
                model.x0(monomials) = [];
                model.variabletype(monomials)=[];
                keptvariables(monomials) = [];
            end
        end
    end
end
