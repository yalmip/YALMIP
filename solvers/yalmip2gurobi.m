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
if any(K.q)
    [F_struc,K,c,Q,UB,LB,x0] = append_normalized_socp(F_struc,K,c,Q,UB,LB,x0);
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

model.params = interfacedata.options.gurobi;
if interfacedata.options.verbose == 0
     model.params.outputflag = 0;
else
     model.params.outputflag = 1;
end

if ~isempty(x0)
    model.start = x0;
end
model.NegativeSemiVar=NegativeSemiVar;