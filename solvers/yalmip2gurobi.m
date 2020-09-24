function model = yalmip2gurobi(interfacedata)

% Retrieve needed data
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

% Recent support for nonconvex case is treated a bit hackish now
% All convex quadratic stuff has been converted to SOCP cones, and then all
% the rest are just kept, with nonlinear monomials remaining in the model
% representation. Save away nonconve info for later, and clean away
nonconvexdata = [];
if any(interfacedata.variabletype) & all(interfacedata.variabletype < 3)
    nonlinearMonoms = find(interfacedata.variabletype);
    s1 = interfacedata.F_struc(1:interfacedata.K.f,1+nonlinearMonoms);
    s2 = interfacedata.F_struc(interfacedata.K.f+1:interfacedata.K.f+interfacedata.K.l,1+nonlinearMonoms);
    k_eq = find(any(s1,2));
    k_ineq = find(any(s2,2));
    if ~isempty(k_eq)
        nonconvexdata.eq = interfacedata.F_struc(k_eq,:);
        interfacedata.F_struc(k_eq,:) = [];
        interfacedata.K.f = interfacedata.K.f - length(k_eq);
    else
        nonconvexdata.eq = [];
    end
    if ~isempty(k_ineq)
        nonconvexdata.ineq = interfacedata.F_struc(interfacedata.K.f + k_ineq,:);
        interfacedata.F_struc(interfacedata.K.f + k_ineq,:) = [];
        interfacedata.K.l = interfacedata.K.l - length(k_ineq);
    else
        nonconvexdata.ineq = [];
    end    
    interfacedata.F_struc(:,1 + nonlinearMonoms) = [];
    interfacedata.c(nonlinearMonoms) = [];
    interfacedata.Q(:,nonlinearMonoms) = [];
    interfacedata.Q(nonlinearMonoms,:) = [];
    interfacedata.lb(nonlinearMonoms) = [];
    interfacedata.ub(nonlinearMonoms) = [];
    if ~isempty(x0)
        x0(nonlinearMonoms) = [];
    end
    F_struc = interfacedata.F_struc;
    c = interfacedata.c;
    Q = interfacedata.Q;    
    lb = interfacedata.lb;
    ub = interfacedata.ub;
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

%[F_struc,K,LB,UB,semicont_variables] = extractSemiContBounds(F_struc,K,LB,UB,semicont_variables);

if ~isempty(semicont_variables)
    % Bounds must be placed in LB/UB
    [LB,UB,cand_rows_eq,cand_rows_lp] = findulb(F_struc,K,LB,UB);
    F_struc(K.f+cand_rows_lp,:)=[];
    F_struc(cand_rows_eq,:)=[];
    K.l = K.l-length(cand_rows_lp);
    K.f = K.f-length(cand_rows_eq);    
    redundant = find(LB<=0 & UB>=0);
    semicont_variables = setdiff(semicont_variables,redundant);        
end

n_original = length(c);
variabletype_original = interfacedata.variabletype;
if any(K.q)
    [F_struc,K,c,Q,UB,LB,x0] = append_normalized_socp(F_struc,K,c,Q,UB,LB,x0);
    if length(c) > interfacedata.variabletype
         interfacedata.variabletype(length(c)) = 0;
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
%[NegativeSemiVar,Q,c,A,lb,ub,semicont_variables] = negateNegativeSemiContVariables(Q,c,A,lb,ub,semicont_variables,[]);
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
            x0(NegativeSemiVar) = -NegativeSemiVar;
        end
    end
end

if ~isempty(K.sos.type)
    for i = 1:length(K.sos.type)
         model.sos(i).index = full(K.sos.variables{i}(:)');
         model.sos(i).weight = full(K.sos.weight{i}(:)');
         if isa(K.sos.type(i),'char')
             model.sos(i).type = str2num(K.sos.type(i));%ns.gurobi.SOS.types(i)= int32(str2num(K.sos.type(i)));
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
        model.quadcon(i).q=zeros(length(c),1);
        model.quadcon(i).rhs=0;
        top = top + n;
    end
end

if ~isempty(nonconvexdata)
    model.params.nonconvex = 2;
    if ~isfield(model,'quadcon')
        model.quadcon = [];
    end
    m = length(model.lb);
    monomials = find(interfacedata.variabletype > 0);
    map = [];
    for j = 1:length(monomials)
        s = find(interfacedata.monomtable(monomials(j),:));
        if length(s) == 1
            map(monomials(j),:) = [s s];
        else
            map(monomials(j),:) = s;
        end
    end
    for i = 1:size(nonconvexdata.eq,1)
        bi = nonconvexdata.eq(i,1);
        row = nonconvexdata.eq(i,2:end);
        qi = row(find(variabletype_original == 0))';
        if length(qi)<m
            % The number of variables has been extended above when SOCPs
            % have been normalized and cast as convex quadratics
            qi(m)=0;
        end
        di = row(monomials);
        Qi = spalloc(m,m,0);
        for k = 1:length(monomials)
            if di(k)
                Qi(map(monomials(k),1),map(monomials(k),2)) = Qi(map(monomials(k),1),map(monomials(k),2)) + di(k)/2;
                Qi(map(monomials(k),2),map(monomials(k),1)) = Qi(map(monomials(k),2),map(monomials(k),1)) + di(k)/2;
            end
        end      
        model.quadcon(end+1).Qc = -Qi;
        model.quadcon(end).q = -qi;
        model.quadcon(end).rhs = bi;
        model.quadcon(end).sense = '=';        
    end
    for i = 1:size(nonconvexdata.ineq,1)
        bi = nonconvexdata.ineq(i,1);
        row = nonconvexdata.ineq(i,2:end);
        qi = row(find(variabletype_original == 0))';
         if length(qi)<m
            % The number of variables has been extended above when SOCPs
            % have been normalized and cast as convex quadratics
            qi(m)=0;
        end
        di = row(monomials);
        Qi = spalloc(m,m,0);
        for k = 1:length(monomials)
            if di(k)
                Qi(map(monomials(k),1),map(monomials(k),2)) = Qi(map(monomials(k),1),map(monomials(k),2)) + di(k)/2;
                Qi(map(monomials(k),2),map(monomials(k),1)) = Qi(map(monomials(k),2),map(monomials(k),1)) + di(k)/2;
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

if isequal(interfacedata.solver.version,'NONCONVEX')
    model.params.nonconvex = 2;
end

if ~isempty(x0)  
    model.start = x0(find(interfacedata.variabletype == 0));
end
model.NegativeSemiVar=NegativeSemiVar;

if isfield(model,'quadcon') && isempty(model.quadcon)
	model = rmfield(model,'quadcon');
end    