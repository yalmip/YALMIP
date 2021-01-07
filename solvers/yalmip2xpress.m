function model = yalmip2xpress(interfacedata)

options = interfacedata.options;
F_struc = interfacedata.F_struc;
H       = interfacedata.Q;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

showprogress('Calling Xpress',options.showprogress);

if ~isempty(semicont_variables)
    % Bounds must be placed in LB/UB
    [lb,ub,cand_rows_eq,cand_rows_lp] = findulb(F_struc,K,lb,ub);
    F_struc(K.f+cand_rows_lp,:)=[];
    F_struc(cand_rows_eq,:)=[];
    K.l = K.l-length(cand_rows_lp);
    K.f = K.f-length(cand_rows_eq);
    redundant = find(lb<=0 & ub>=0);
    semicont_variables = setdiff(semicont_variables,redundant);
end

% Notation used
f = c;
A = -F_struc(:,2:end);
b = F_struc(:,1);
rtype = repmat('L',size(F_struc,1),1);
rtype(1:K.f) = 'E';

% XPRESS assumes semi-continuous variables only can take positive values so
% we negate semi-continuous violating this
NegativeSemiVar = [];
if ~isempty(semicont_variables)
    NegativeSemiVar = find(lb(semicont_variables) < 0);
    if ~isempty(NegativeSemiVar)
        temp = ub(semicont_variables(NegativeSemiVar));
        ub(semicont_variables(NegativeSemiVar)) = -lb(semicont_variables(NegativeSemiVar));
        lb(semicont_variables(NegativeSemiVar)) = -temp;
        if ~isempty(A)
            A(:,semicont_variables(NegativeSemiVar)) = -A(:,semicont_variables(NegativeSemiVar));
        end
        f(semicont_variables(NegativeSemiVar)) = -f(semicont_variables(NegativeSemiVar));
        H(:,semicont_variables(NegativeSemiVar)) = -H(:,semicont_variables(NegativeSemiVar));
        H(semicont_variables(NegativeSemiVar),:) = -H(semicont_variables(NegativeSemiVar),:);
    end
end
clim = zeros(length(f),1);
clim(semicont_variables) = lb(semicont_variables);

ctype = char(ones(length(f),1)*67);
ctype(setdiff(integer_variables,semicont_variables)) = 'I';
ctype(binary_variables)  = 'B';
ctype(setdiff(semicont_variables,integer_variables)) = 'S';
ctype(intersect(semicont_variables,integer_variables)) = 'R';

options.xpress = setupOptions(options);
if options.verbose == 0
    options.xpress.OUTPUTLOG = 0;
    options.xpress.MIPLOG = 0;
    options.xpress.LPLOG = 0;
end

sos = [];
if ~isempty(K.sos.type)
    for i = 1:length(K.sos.weight)
        sos(i).type = K.sos.type(i);
        sos(i).ind =  K.sos.variables{i}-1;
        sos(i).wt = (1:length(K.sos.weight{i}))';
    end
end

if size(A,1)==0
    A = zeros(1,length(f));
    b = 1;
end

model.H = 2*H;
model.f = f;
model.A = A;
model.b = b;
model.lb = lb;
model.ub = ub;
model.ctype = ctype;
model.rtype = rtype;
model.clim = clim;
model.sos = sos;
model.ops = options.xpress;
model.extra.semicont_variables = semicont_variables;
model.extra.integer_variables = integer_variables;
model.extra.binary_variables = binary_variables;
model.extra.NegatedSemiVar = NegativeSemiVar;

function ops = setupOptions(options);

ops = options;
if options.verbose == 0
    ops.OUTPUTLOG = 0;
    ops.MIPLOG = 0;
    ops.LPLOG = 0;
end

% cNames = fieldnames(options.xpress);
% for i = 1:length(cNames)
%     s = getfield(options.xpress,cNames{i});
%     if ~isempty(s)
%         ops = setfield(ops,cNames{i},s);
%     end
% end