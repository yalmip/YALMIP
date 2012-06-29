function model = yalmip2gurobimex(interfacedata)
% Retrieve needed data
options = interfacedata.options;
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
interfacedata.gettime = 0;
n = length(c);

if ~isempty(ub)
    LB = lb;
    UB = ub;
    LB(binary_variables)  = round(LB(binary_variables));
    UB(binary_variables)  = round(UB(binary_variables));
    LB(integer_variables) = round(LB(integer_variables));
    UB(integer_variables) = round(UB(integer_variables));
else
    LB = [];
    UB = [];
end

if options.showprogress;showprogress('Calling GUROBI',options.showprogress);end

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

SENSE = 1;     % Minimize
C = full(c);   % Must be full
if size(F_struc,1)>0
    B = full(F_struc(:,1));         % Must be full
    A =-F_struc(:,2:end);
else
    B = [];
    A = [];
end

% Optimized code, make a lot of difference when you make this call 10000
% times in a branch and bound setting...
CTYPE = [char(ones(K.f,1)*61); char(ones(K.l,1)*60)];
VARTYPE = char(ones(length(c),1)*67);
VARTYPE(setdiff(integer_variables,semicont_variables)) = 'I';
VARTYPE(binary_variables)  = 'B';  % Should not happen except from bmibnb
VARTYPE(setdiff(semicont_variables,integer_variables)) = 'S';  % Should not happen except from bmibnb
VARTYPE(intersect(semicont_variables,integer_variables)) = 'N';

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
        C(semicont_variables(NegativeSemiVar)) = -C(semicont_variables(NegativeSemiVar));
        if ~isempty(x0)
            x0(NegativeSemiVar) = -NegativeSemiVar;
        end
    end
end

if nnz(Q)>0
    [ii,jj,kk] = find(Q);
    ii = ii-1;
    jj = jj-1;
    comp = computer;
    % According to Wotao's testing
    if isempty(strfind(comp,'64')) | isequal(comp,'GLNXA64') | isequal(comp,'PCWIN64')
        options.gurobi.QP.qrow = int32(ii)';
        options.gurobi.QP.qcol = int32(jj)';
    elseif strfind(comp,'MACI64')
        options.gurobi.QP.qrow = int32(ii)';
        options.gurobi.QP.qcol = int32(jj)';    
    else
        options.gurobi.QP.qrow = int64(ii)';
        options.gurobi.QP.qcol = int64(jj)';
    end
    options.gurobi.QP.qval = kk';
end
        
if ~options.verbose
    options.gurobi.DisplayInterval = 0;
    options.gurobi.Display = 0;
end

if ~isempty(x0)
    options.gurobi.Start = x0(:)';
end

if options.savedebug
    save gurobidebug
end

if ~isempty(K.sos.type)
    options.gurobi.SOS.weights = spalloc(length(c),length(K.sos.type),0);
    for i = 1:length(K.sos.type)
        options.gurobi.SOS.types(i)= int32(str2num(K.sos.type(i)));
        options.gurobi.SOS.weights(K.sos.variables{i},i) = K.sos.weight{i};
    end
end
model.C = C;