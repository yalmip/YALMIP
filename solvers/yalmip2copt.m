function model = yalmip2copt(interfacedata)
F_struc            = interfacedata.F_struc;
K                  = interfacedata.K;
c                  = interfacedata.c;
lb                 = interfacedata.lb;
ub                 = interfacedata.ub;
x0                 = interfacedata.x0;
integer_variables  = interfacedata.integer_variables;
binary_variables   = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;
n                  = length(c);

nlindim = K.f + K.l;
if ~isempty(K.q) && any(K.q)
    nsocdim = sum(K.q); % Standard quadratic cone only
else
    nsocdim = 0;
end
if ~isempty(K.e) && K.e > 0
    nexpdim = 3 * K.e;
else
    nexpdim = 0;
end

if ~isempty(K.s) && any(K.s)
    data.A = -F_struc(:, 2:end);
    data.b = full(F_struc(:, 1));
    data.c = full(c);
    
    dims = [];
    dims.f = K.f;
    dims.l = K.l;
    if ~isempty(K.q)
        dims.q = K.q;
    end
    if ~isempty(K.e)
        dims.ep = K.e;
    end
    dims.s = K.s;
    
    if nexpdim > 0
        expA = data.A(1 + nlindim + nsocdim:nlindim + nsocdim + nexpdim, :);
        expb = data.b(1 + nlindim + nsocdim:nlindim + nsocdim + nexpdim, :);
    end
    sdpA = data.A(1 + nlindim + nsocdim + nexpdim:end, :);
    sdpb = data.b(1 + nlindim + nsocdim + nexpdim:end, :);
    
    data.A = data.A(1:nlindim + nsocdim, :);
    data.b = data.b(1:nlindim + nsocdim, :);
    
    top = 1;
    if nexpdim > 0
        for i = 1:K.e
            data.A = [data.A; expA(top + [2; 1; 0], :)];
            data.b = [data.b; expb(top + [2; 1; 0], :)];
            top = top + 3;
        end
    end
    
    top = 1;
    for i = 1:length(K.s)
        A = sdpA(top:top + K.s(i)^2 - 1, :);
        b = sdpb(top:top + K.s(i)^2 - 1, :);
        n = K.s(i);
        ind = find(tril(ones(n)));
        A = A(ind, :);
        b = b(ind);
        data.A = [data.A; A];
        data.b = [data.b; b];
        top = top + K.s(i)^2;
    end
    
    conedata.A = data.A';
    conedata.c = -data.b;
    conedata.b = -data.c;
    conedata.K = dims;
    conedata.objsen = 'max';
    model.conedata = conedata;
    
    model.params = interfacedata.options.copt;
    if interfacedata.options.verbose == 0
        model.params.Logging = 0;
    else
        model.params.Logging = 1;
    end
    
    return
end

if ~isempty(ub)
    LB = lb;
    UB = ub;
    if ~isempty(binary_variables)
        LB(binary_variables) = round(LB(binary_variables));
        UB(binary_variables) = round(UB(binary_variables));
    end
    if ~isempty(integer_variables)
        LB(integer_variables) = round(LB(integer_variables));
        UB(integer_variables) = round(UB(integer_variables));
    end
else
    LB = -inf(n, 1);
    UB = +inf(n, 1);
end

if size(F_struc, 1) > 0
    A   = -F_struc(:, 2:end);
    RHS = full(F_struc(:, 1));
    LHS = -inf(length(RHS), 1);
else
    A   = sparse(ones(1, length(c)));
    RHS = +inf(length(c), 1);
    LHS = -inf(length(c), 1);
end

if K.f > 0
    LHS(1:K.f) = RHS(1:K.f);
end

VARTYPE = char(ones(length(c), 1) * 67);
if isempty(semicont_variables)
    VARTYPE(integer_variables) = 'I';
    VARTYPE(binary_variables)  = 'B';
end

model.objcon = full(interfacedata.f);
if isempty(A)
    model.A = spalloc(0, length(c), 0);
else
    model.A = sparse(A);
end
model.obj   = full(c);
if any(interfacedata.Q)
    model.Q     = interfacedata.Q;
end
model.lb    = LB;
model.ub    = UB;
model.vtype = VARTYPE;
model.lhs   = LHS;
model.rhs   = RHS;

if ~isempty(x0)
    model.start = x0;
end

norigcol = size(F_struc, 2) - 1;
if nsocdim > 0
    if nlindim > 0
        model.A = [model.A, [spalloc(nlindim, nsocdim, 0); speye(nsocdim)]];
    else
        model.A = [model.A, speye(nsocdim)];
    end
    model.obj   = [model.obj; zeros(nsocdim, 1)];
    model.lb    = [model.lb; -inf(nsocdim, 1)];
    model.ub    = [model.ub; +inf(nsocdim, 1)];
    model.vtype = [model.vtype; char(ones(nsocdim, 1) * 67)];
    
    top = norigcol;
    for i = 1:length(K.q)
        model.cone(i).type = 1;
        model.cone(i).vars = top + 1:top + K.q(i);
        top = top + K.q(i);
    end
end

if nexpdim > 0
    expmask = zeros(nexpdim, 1);
    for i = 1:K.e
        expmask(3*i - 2) = 3*i;
        expmask(3*i - 1) = 3*i - 1;
        expmask(3*i)     = 3*i - 2;
    end
    eyemat = speye(nexpdim);
    expmat = eyemat(expmask, :);
    
    if nlindim + nsocdim > 0
        model.A = [model.A, [spalloc(nlindim + nsocdim, nexpdim, 0); expmat]];
    else
        model.A = [model.A, expmat];
    end
    model.obj   = [model.obj; zeros(nexpdim, 1)];
    model.lb    = [model.lb; -inf(nexpdim, 1)];
    model.ub    = [model.ub; +inf(nexpdim, 1)];
    model.vtype = [model.vtype; char(ones(nexpdim, 1) * 67)];

    top = norigcol + nsocdim;
    for i = 1:K.e
        model.expcone(i).type = 3;
        model.expcone(i).vars = [top + 1; top + 2; top + 3];
        top = top + 3;
    end
end

if nsocdim > 0 || nexpdim > 0
    if any(interfacedata.Q)
        model.Q(length(model.obj),length(model.obj)) = 0;
    end
    if ~isempty(x0)
        model.start = [model.start; zeros(nsocdim + nexpdim, 1)];
    end
    model.lhs(1 + nlindim:end) = model.rhs(1 + nlindim:end);
end

if ~isempty(K.sos.type)
    for i = 1:length(K.sos.type)
        if isa(K.sos.type(i), 'char')
            model.sos(i).type = str2num(K.sos.type(i));
        else
            model.sos(i).type = K.sos.type(i);
        end
        model.sos(i).vars    = full(K.sos.variables{i}(:)');
        model.sos(i).weights = full(K.sos.weight{i}(:)');
    end
end

model.params = interfacedata.options.copt;
if interfacedata.options.verbose == 0
    model.params.Logging = 0;
else
    model.params.LogFile = 'copt.log';
end

