function model = yalmip2gurobinonlinear(interfacedata)

% Experimental support for nonlinear stuff
% Simplifies by keeping every monomial (so lots of zeros in final model)
F_struc = interfacedata.F_struc;
Q       = interfacedata.Q;
c       = interfacedata.c;
K       = interfacedata.K;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
x0      = interfacedata.x0;
n = length(c);

% Model might have been compiled for the nonlinear layer which means we
% have to compile the Q-part of the objective here instead
quadraticMonoms = find(interfacedata.variabletype <= 2 & interfacedata.variabletype >= 1);
if ~isempty(quadraticMonoms)
    for j = quadraticMonoms(find(interfacedata.c(quadraticMonoms)))
        cj = interfacedata.c(j);
        interfacedata.c(j) = 0;        
        i = find(interfacedata.monomtable(j,:));
        if interfacedata.variabletype(j)==1        
            interfacedata.Q(i(1),i(2)) = cj/2;
            interfacedata.Q(i(2),i(1)) = cj/2;
        else            
            interfacedata.Q(i,i) = cj;
        end        
    end
    c = interfacedata.c;
    Q = interfacedata.Q;
end

if ~isempty(interfacedata.evalMap) && any(interfacedata.variabletype)
    % We might have the case such as exp(x^2). YALMIP handles this natively
    % by computing composite stuff for nonlinear solvers when argument is a
    % simple monomial. Gurobi on the other hand can only have linear variab
    % les. Hence, we have to introduce a new variable
    for i = 1:length(interfacedata.evalMap)
        if interfacedata.variabletype(interfacedata.evalMap{i}.variableIndex)
            % Introduce a new linear variable
            n = size(interfacedata.monomtable,1)+1;
            interfacedata.monomtable(n,n) = 1;
            interfacedata.variabletype(n) = 0;
            % Use as argument
            old = interfacedata.evalMap{i}.variableIndex;
            interfacedata.evalMap{i}.variableIndex = n;
            interfacedata.c(n) = 0;c = interfacedata.c;
            interfacedata.Q(n,n) = 0;Q = interfacedata.Q;
            if ~isempty(interfacedata.F_struc);interfacedata.F_struc(1,n+1) = 0;end
            if ~isempty(lb);lb(n)=-inf;end
            if ~isempty(ub);ub(n)=inf;end
            row = zeros(1,n);row(n)=1;row(old) = -1;
            interfacedata.F_struc = [0 row;interfacedata.F_struc];
            interfacedata.K.f = interfacedata.K.f + 1;
        end
    end
end

% All convex quadratics have been converted to SOCP cones, and then all
% the rest involving quadratics are just kept. Find these and put in
% gurobis dedicated format
quadraticconstraints = [];
%quadraticMonoms = find(interfacedata.variabletype <= 2 & interfacedata.variabletype >= 1);
if ~isempty(quadraticMonoms)
    equalityRows = interfacedata.F_struc(1:interfacedata.K.f,1+quadraticMonoms);
    inequalityRows = interfacedata.F_struc(interfacedata.K.f+1:interfacedata.K.f+interfacedata.K.l,1+quadraticMonoms);
    k_eq = find(any(equalityRows,2));
    k_ineq = find(any(inequalityRows,2));
    if ~isempty(k_eq)
        quadraticconstraints.eq = interfacedata.F_struc(k_eq,:);
        interfacedata.F_struc(k_eq,:) = [];
        interfacedata.K.f = interfacedata.K.f - length(k_eq);
    else
        quadraticconstraints.eq = [];
    end
    if ~isempty(k_ineq)
        quadraticconstraints.ineq = interfacedata.F_struc(interfacedata.K.f + k_ineq,:);
        interfacedata.F_struc(interfacedata.K.f + k_ineq,:) = [];
        interfacedata.K.l = interfacedata.K.l - length(k_ineq);
    else        
        quadraticconstraints.ineq = [];
    end    
    if isempty(k_eq) && isempty(k_ineq)
        quadraticconstraints = [];
    end
    F_struc = interfacedata.F_struc;
    K = interfacedata.K;
end

if ~isempty(ub)
    LB = lb;
    UB = ub;
    if ~isempty(binary_variables)
        LB(binary_variables)  = round(LB(binary_variables));
        UB(binary_variables)  = round(UB(binary_variables));
    end
    if ~isempty(integer_variables)
        LB(integer_variables) = round(LB(integer_variables));
        UB(integer_variables) = round(UB(integer_variables));
    end
else
    LB = -inf(n,1);
    UB = inf(n,1);
end

if ~isempty(semicont_variables)
    % Bounds must be placed in LB/UB
    [LB,UB,cand_rows_eq,cand_rows_lp] = find_lp_bounds(F_struc,K,LB,UB);
    F_struc(K.f+cand_rows_lp,:)=[];
    F_struc(cand_rows_eq,:)=[];
    K.l = K.l-length(cand_rows_lp);
    K.f = K.f-length(cand_rows_eq);
    redundant = find(LB<=0 & UB>=0);
    semicont_variables = setdiff(semicont_variables,redundant);
end

n_original = length(c);
if any(K.q)
    [F_struc,K,c,Q,UB,LB,x0] = append_normalized_socp(F_struc,K,c,Q,UB,LB,x0);
    if length(c) > interfacedata.variabletype
        interfacedata.variabletype(length(c)) = 0;
        interfacedata.monomtable(end,length(c))=0;
    end
end

if size(F_struc,1)>0
    B = full(F_struc(:,1));         % Must be full
    A =-F_struc;
    A(:,1)=[];
else
    B = [];
    A = [];
end

% Optimized code, make a lot of difference when you make this call 10000
% times in a branch and bound setting...
CTYPE = char([ones(K.f,1)*61; ones(K.l,1)*60]);
VARTYPE = char(ones(length(c),1)*67);
if isempty(semicont_variables)
    VARTYPE(integer_variables) = 'I';
    VARTYPE(binary_variables)  = 'B';  % Should not happen except from bmibnb
else
    VARTYPE(setdiff(integer_variables,semicont_variables)) = 'I';
    VARTYPE(binary_variables)  = 'B';  % Should not happen except from bmibnb
    VARTYPE(setdiff(semicont_variables,integer_variables)) = 'S';  % Should not happen except from bmibnb
    VARTYPE(intersect(semicont_variables,integer_variables)) = 'N';
end

% Gurobi assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
NegativeSemiVar = [];
if ~isempty(semicont_variables)
    NegativeSemiVar = find(UB(semicont_variables) < 0);
    if ~isempty(NegativeSemiVar)
        temp = UB(semicont_variables(NegativeSemiVar));
        UB(semicont_variables(NegativeSemiVar)) = -LB(semicont_variables(NegativeSemiVar));
        LB(semicont_variables(NegativeSemiVar)) = -temp;
        A(:,semicont_variables(NegativeSemiVar)) = -A(:,semicont_variables(NegativeSemiVar));
        c(semicont_variables(NegativeSemiVar)) = -c(semicont_variables(NegativeSemiVar));
        if ~isempty(x0)
            x0(semicont_variables(NegativeSemiVar)) = -x0(semicont_variables(NegativeSemiVar));
        end
    end
end

if ~isempty(K.sos.type)
    for i = 1:length(K.sos.type)
        model.sos(i).index = full(K.sos.variables{i}(:)');
        model.sos(i).weight = full(K.sos.weight{i}(:)');
        if isa(K.sos.type(i),'char')
            model.sos(i).type = str2num(K.sos.type(i));
        else
            model.sos(i).type = K.sos.type(i);
        end
    end
end

if isempty(A)
    model.A = spalloc(0,length(c),0);
else
    model.A = sparse(A);
end
model.obj = full(c);
model.sense = CTYPE;
model.rhs = full(B);
model.lb = LB;
model.ub = UB;
model.objcon = full(interfacedata.f);
model.vtype = VARTYPE;
model.Q = sparse(Q);
model.params = interfacedata.options.gurobi;

if ~isequal(K.q,0)
    top = n_original + 1;
    for i = 1:length(K.q)
        n = K.q(i);
        Qi = sparse(top:top+n-1,top:top+n-1,[-1 repmat(1,1,n-1)],length(c),length(c));
        model.quadcon(i).Qc=Qi;
        model.quadcon(i).q=sparse(length(c),1);
        model.quadcon(i).rhs=0;
        top = top + n;
    end
end

if ~isempty(quadraticconstraints)
    model.params.nonconvex = 2;
    if ~isfield(model,'quadcon')
        model.quadcon = [];
    end
    m = length(model.lb);    
              
    map = [];
    for j = 1:length(quadraticMonoms)
        s = find(interfacedata.monomtable(quadraticMonoms(j),:));
        if length(s) == 1
            map(quadraticMonoms(j),:) = [s s];
        else
            map(quadraticMonoms(j),:) = s;
        end
    end
    
    for i = 1:size(quadraticconstraints.eq,1)
        bi = quadraticconstraints.eq(i,1);
        row = quadraticconstraints.eq(i,2:end);
        di = row(quadraticMonoms);
        qi = row(:);qi(quadraticMonoms) = 0;        
        if length(qi)<m
            % The number of variables has been extended above when SOCPs
            % have been normalized and cast as convex quadratics
            qi(m)=0;
        end
        Qi = spalloc(m,m,0);
        for k = find(di)
            Qi(map(quadraticMonoms(k),1),map(quadraticMonoms(k),2)) = Qi(map(quadraticMonoms(k),1),map(quadraticMonoms(k),2)) + di(k)/2;
            Qi(map(quadraticMonoms(k),2),map(quadraticMonoms(k),1)) = Qi(map(quadraticMonoms(k),2),map(quadraticMonoms(k),1)) + di(k)/2;
        end
        model.quadcon(end+1).Qc = -Qi;
        model.quadcon(end).q = -qi;
        model.quadcon(end).rhs = bi;
        model.quadcon(end).sense = '=';
    end
    for i = 1:size(quadraticconstraints.ineq,1)
        bi = quadraticconstraints.ineq(i,1);
        row = quadraticconstraints.ineq(i,2:end);
        di = row(quadraticMonoms);
        qi = row(:);qi(quadraticMonoms) = 0;          
        if length(qi)<m
            % The number of variables has been extended above when SOCPs
            % have been normalized and cast as convex quadratics
            qi(m)=0;
        end        
        Qi = spalloc(m,m,0);
        for k = 1:length(quadraticMonoms)
            if di(k)
                Qi(map(quadraticMonoms(k),1),map(quadraticMonoms(k),2)) = Qi(map(quadraticMonoms(k),1),map(quadraticMonoms(k),2)) + di(k)/2;
                Qi(map(quadraticMonoms(k),2),map(quadraticMonoms(k),1)) = Qi(map(quadraticMonoms(k),2),map(quadraticMonoms(k),1)) + di(k)/2;
            end
        end
        model.quadcon(end+1).Qc = -Qi;
        model.quadcon(end).q = -qi;
        model.quadcon(end).rhs = bi;
        model.quadcon(end).sense = '<';
    end
end

if interfacedata.options.verbose == 0
    
    model.params.outputflag = 0;
else
    model.params.outputflag = 1;
end

if isequal(interfacedata.solver.version,'NONCONVEX') || isequal(interfacedata.solver.version,'NONLINEAR')
    model.params.nonconvex = 2;
end

% YALMIP defaults to outer approximation
if ~isfield(interfacedata.options.gurobi,'FuncNonLinear')
    model.params.FuncNonLinear = 1;
end

if ~isempty(x0)
    model.start = x0;
    if length(model.start)~=length(model.obj)
        model.start=[];
    end
end
model.NegativeSemiVar=NegativeSemiVar;

if isfield(model,'quadcon') && isempty(model.quadcon)
    model = rmfield(model,'quadcon');
end

% Monomials of higher order
high_monoms = find(interfacedata.variabletype > 2);
if ~isempty(high_monoms)
    model.genconpow = [];
    for i = high_monoms(:)'
        monoms = interfacedata.monomtable(i,:);
        [~,var,p] = find(monoms);
        if length(p)==1
            % Variable x(i) is x(var)^p
            model.genconpow = [model.genconpow struct('xvar',var,'yvar', i,'a',p)];
        elseif length(p)>1 && interfacedata.variabletype(i)>2
            % Ocuh, monomial term with multiple variables not supported
            error('GUROBI does not support general polynomial expressions')
        end
    end
end

% Nonlinear operators
if length(interfacedata.evalMap) > 0
    
    functions = {'exp','log','sin','cos','tan','logistic'};
    for i = 1:length(functions)
        if any(cellfun(@(x) strcmp(x.fcn, functions{i}), interfacedata.evalMap))
            eval(['model.gencon' functions{i} '=[];']);
        end        
    end
    if any(cellfun(@(x) strcmp(x.fcn, 'power_internal1'), interfacedata.evalMap))
        model.genconexpa = [];
    end
        
    for i = 1:length(interfacedata.evalMap)
        f = interfacedata.evalMap{i};
        switch f.fcn
            case 'log'
                model.genconlog = [model.genconlog struct('xvar',f.variableIndex,'yvar', f.computes)];
            case 'exp'
                model.genconexp = [model.genconexp struct('xvar',f.variableIndex,'yvar', f.computes)];
            case 'sin'
                model.genconsin = [model.genconsin struct('xvar',f.variableIndex,'yvar', f.computes)];
            case 'cos'
                model.genconcos = [model.genconcos struct('xvar',f.variableIndex,'yvar', f.computes)];
            case 'tan'
                model.gencontan = [model.gencontan struct('xvar',f.variableIndex,'yvar', f.computes)];                
            case 'logistic'
                model.genconlogistic = [model.genconlogistic struct('xvar',f.variableIndex,'yvar', f.computes)];                
            case 'power_internal1'
                % a^x
                aux = struct('xvar',f.variableIndex,'yvar', f.computes,'a',interfacedata.evalMap{i}.arg{2});
                model.genconexpa = [model.genconexpa aux];                
            otherwise
                error(['GUROBI does not support "'  f.fcn '"']);
        end
    end
end